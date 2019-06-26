"""
Defines :class:`StructUnit` classes.

:class:`StructUnit` represents the monomers that make up
macromolecules. These are commonly refered to as ``building blocks`` in
the documentation. :class:`StructUnit` holds information concerning
only a single building block molecule. For example, the number of atoms
and bonds a building block may have. It also has information about the
functional groups present in the building block molecule (see
:class:`.FGInfo`). The class also allows manipulation of the lone
building block molecule, such as rotations and translations.

The :class:`StructUnit` should be inherited as necessary. For example,
:class:`StructUnit2` adds manipulations relavant to molecules with 2
functional groups. :class:`StructUnit3` adds manipulations relavant to
molecules with 3 or more functional groups. If you have a monomer which
needs specific information or manipulations, give it its own class.

:class:`.StructUnit` contains :class:`.FunctionalGroup` instances
representing functional groups used during macromolecular assembly.
There are divided into bonders and deleters (see the documentation of
:class:`.FGInfo` and :class:`.FunctionalGroup`), which determines
which atoms form bonds and which are removed during assembly.

"""


import logging
import os
import numpy as np
import itertools as it
import rdkit.Chem.AllChem as rdkit
from rdkit import DataStructs
from glob import glob
from functools import partial
from scipy.spatial.distance import euclidean
from collections import defaultdict
from inspect import signature

from .molecule import Molecule
from ..functional_groups import FunctionalGroup
from ..functional_groups import functional_group_infos as fg_infos
from ..functional_groups import functional_groups as fgs
from ...utilities import (flatten,
                          normalize_vector,
                          rotation_matrix,
                          vector_theta,
                          mol_from_mae_file,
                          rotation_matrix_arbitrary_axis,
                          remake,
                          dedupe,
                          OPTIONS)


logger = logging.getLogger(__name__)


class CachedStructUnit(type):
    """
    A metaclass for making :class:`StructUnit` create cached instances.

    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.cache = dict()

    def __call__(self, *args, **kwargs):
        # Get the arguments given to the initializer as a dictionary
        # mapping argument name to argument value.
        sig = signature(self.__init__)
        sig = sig.bind_partial(self, *args, **kwargs)
        sig.apply_defaults()
        sig = sig.arguments

        is_str = isinstance(sig['mol'], str)
        if is_str:
            if os.path.exists(sig['mol']):
                _, ext = os.path.splitext(sig['mol'])

                if ext not in self.init_funcs:
                    raise TypeError(
                        f'Unable to initialize from "{ext}" files.'
                    )
                mol = remake(self.init_funcs[ext](sig['mol']))

            else:
                mol = remake(
                    rdkit.MolFromMolBlock(molBlock=sig['mol'],
                                          sanitize=False,
                                          removeHs=False)
                )

        elif isinstance(sig['mol'], rdkit.Mol):
            mol = remake(sig['mol'])

        # Get the name of the functional groups provided to the
        # initializer or get it from the path.
        if sig['functional_groups'] is not None:
            functional_groups = sig['functional_groups']

        elif is_str and os.path.exists(sig['mol']):
            functional_groups = tuple((
                fg.name for fg in fgs if fg.name in sig['mol']
            ))

        else:
            functional_groups = ()

        key = self.gen_key(mol, functional_groups)
        if key in self.cache and OPTIONS['cache']:
            return self.cache[key]
        else:
            obj = super().__call__(*args, **kwargs)
            obj.key = key
            if OPTIONS['cache']:
                self.cache[key] = obj
            return obj


class StructUnit(Molecule, metaclass=CachedStructUnit):
    """
    Represents the building blocks of macromolecules.

    The goal of this class is to conveniently store information about,
    and perform operations on, single instances of macromolecular
    building blocks.

    Attributes
    ----------
    init_funcs : :class:`dict`
        This dictionary holds the various functions which can be used
        to initialize ``rdkit`` molecules and pairs them with the
        appropriate file extension.

    file : :class:`str`
        The full path to the molecular structure file holding the
        molecule. The supported file formats are the keys in
        :attr:`init_funcs`. As long as a file with one of these
        extensions is provided, the correct initialization function
        will be used.

    func_groups : :class:`tuple` of :class:`.FunctionalGroup`
        The functional groups present in the molecule. The id of
        each :class:`.FunctionalGroup` should match its index.

    func_group_infos : :class:`tuple` of :class:`.FGInfo`
        The functional group types present in the molecule.

    """

    init_funcs = {'.mol': partial(rdkit.MolFromMolFile,
                                  sanitize=False, removeHs=False),

                  '.sdf': partial(rdkit.MolFromMolFile,
                                  sanitize=False, removeHs=False),

                  '.mol2': partial(rdkit.MolFromMol2File,
                                   sanitize=False, removeHs=False),

                  '.mae': mol_from_mae_file,

                  '.pdb': partial(rdkit.MolFromPDBFile,
                                  sanitize=False, removeHs=False)}

    def __init__(self, mol, functional_groups=None, name="", note=""):
        """
        Initializes a :class:`StructUnit` instance.

        Parameters
        ----------
        mol : :class:`str` or :class:`rdkit.Mol`
            Can be one of 3 things:

                1. A path to a molecular structure file.
                2. A :class:`rdkit.Mol` object.
                3. V3000 MDL Mol block.

        functional_groups : :class:`list` of :class:`str`, optional
            The names of the functional groups which are to have atoms
            tagged. If ``None``, a functional group name found in the
            path `file`  is used. If no functional groups are provided
            to this parameter and the name of one is not present in
            `file`, no tagging is done.

        name : :class:`str`, optional
            A name which can be optionally given to the molecule for
            easy identification.

        note : :class:`str`, optional
            A note or comment about the molecule.

        """

        self.file = None

        if isinstance(mol, str):
            if os.path.exists(mol):
                self.file = mol
                _, ext = os.path.splitext(mol)

                if ext not in self.init_funcs:
                    raise TypeError(
                        f'Unable to initialize from "{ext}" files.'
                    )
                self.mol = remake(self.init_funcs[ext](mol))

            else:
                self.mol = remake(
                    rdkit.MolFromMolBlock(molBlock=mol,
                                          removeHs=False,
                                          sanitize=False)
                )

        elif isinstance(mol, rdkit.Mol):
            self.mol = remake(mol)

        # Update the property cache of each atom. This updates things
        # like valence.
        for atom in self.mol.GetAtoms():
            atom.UpdatePropertyCache()

        # If no functional group names passed, check if any functional
        # group names appear in the file path.
        if functional_groups is None:
            if isinstance(mol, str) and os.path.exists(mol):
                functional_groups = [
                    fg.name for fg in fgs if fg.name in mol
                ]
            else:
                functional_groups = []

        self.func_groups = tuple(
            self.functional_groups(functional_groups)
        )

        self.func_group_infos = tuple(dedupe(
            iterable=(fg.info for fg in self.func_groups),
            key=lambda info: info.name
        ))

        super().__init__(name, note)

    @classmethod
    def init_random(cls, db, functional_groups=None, name="", note=""):
        """
        Picks a random file from `db` to initialize from.

        Parameters
        ----------
        db : :class:`str`
            A path to a database of molecular files.

        functional_groups : :class`list` of :class:`str`, optional
            The name of a functional groups which the molecules in `db`
            have. By default it is assumed the name is present in the
            path of the files.

        name : :class:`str`, optional
            The name to be given to the created molecule.

        note : :class:`str`, optional
            A note to be given to the created molecule.

        Returns
        -------
        :class:`StructUnit`
            A random molecule from `db`.

        Raises
        ------
        :class:`RuntimeError`
            If no files in `db` could be initialized from.

        """

        files = glob(os.path.join(db, '*'))
        np.random.shuffle(files)

        for molfile in files:
            try:
                return cls(molfile, functional_groups, name, note)

            except Exception:
                msg = (f'Could not initialize '
                       f'{cls.__name__} from {molfile}.')
                logger.warning(msg)
        raise RuntimeError(
            f'No files in "{db}" could be initialized from.'
        )

    def bonder_centroid(self, conformer=-1):
        """
        Returns the centroid of the bonder atoms.

        Parameters
        ----------
        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        :class:`numpy.ndarray`
            An array holding the midpoint of the bonder atoms.

        """

        # Take the centroid of the bonder centroids.
        nfgs = len(self.func_groups)
        return sum(self.bonder_centroids(conformer)) / nfgs

    def bonder_centroids(self, conformer=-1):
        """
        Calculates the centriod of bonder atoms in each fg.

        The centroids are yielded in order, with the centroid of
        the functional group with an id of ``0`` first.

        Parameters
        ----------
        conformer : :class:`int`, optional
            The conformer to use.

        Yields
        ------
        :class:`numpy.ndarray`
            The bonder centroid of a functional group.

        """

        for fg in self.func_groups:
            bonder_coords = (self.atom_coords(bonder, conformer)
                             for bonder in fg.bonder_ids)
            yield sum(bonder_coords) / len(fg.bonder_ids)

    def all_bonder_distances(self, conformer=-1):
        """
        Yield distances between all pairs of bonder centroids.

        Parameters
        ----------
        conformer : :class:`int`, optional
            The conformer to use.

        Yields
        ------
        :class:`tuple`
            A :class:`tuple` of the form ``(3, 54, 12.54)``.
            The first two elements are the ids of the involved
            functional groups and the third element is the distance
            between them.

        """

        centroids = it.combinations(
                        iterable=self.bonder_centroids(conformer),
                        r=2)
        ids = it.combinations(iterable=range(len(self.func_groups)),
                              r=2)
        for (id1, id2), (c1, c2) in zip(ids, centroids):
            yield id1, id2, euclidean(c1, c2)

    def bonder_direction_vectors(self, conformer=-1):
        """
        Yields the direction vectors between bonder atoms.

        First the centroids of bonder atoms within the same functional
        group are found. Then direction vectors between all pairs of
        these centroids are yielded. Each pair is only yielded once.

        The yielded vector is normalized.

        Parameters
        ----------
        conformer : :class:`int`, optional
            The conformer to use.

        Yields
        ------
        :class:`tuple`
            They yielded tuple has the form

            .. code-block:: python

                (3, 54, np.array([12.2, 43.3, 9.78]))

            The first two elements of the tuple represent the ids
            of the start and end fgs of the vector, respectively. The
            array is the direction vector running between the
            functional group positions.

        """

        centroids = it.combinations(
                        iterable=self.bonder_centroids(conformer),
                        r=2)
        ids = it.combinations(
                  iterable=range(len(self.func_groups)),
                  r=2)
        for (id1, id2), (c1, c2) in zip(ids, centroids):
            yield id2, id1, normalize_vector(c1-c2)

    def bonder_position_matrix(self, conformer=-1):
        """
        Returns a matrix holding the positions of bonder centroids.

        Parameters
        ----------
        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        :class:`numpy.ndarray`
            The array has the shape ``[3, n]``. Each column holds the
            x, y and z coordinates of a bonder centroid. The index of
            the column corresponds to the id of the functional group.

        """

        return np.array(list(self.bonder_centroids(conformer))).T

    def centroid_centroid_dir_vector(self, conformer=-1):
        """
        Returns the direction vector between the 2 molecular centroids.

        The first molecular centroid is the centroid of the entire
        molecule. The second molecular centroid is given by
        :meth:`bonder_centroid`.

        Parameters
        ---------
        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        :class:`numpy.ndarray`
            The normalized direction vector running from the centroid
            of the bonder atoms to the molecular centroid.

        """

        # If the bonder centroid and centroid are in the same position,
        # the centroid - centroid vector should be orthogonal to the
        # bonder direction vector.
        close = np.allclose(a=self.centroid(conformer),
                            b=self.bonder_centroid(conformer),
                            atol=1e-5)
        if close:
            *_, bvec = next(self.bonder_direction_vectors(conformer))
            # Construct a secondary vector by finding the minimum
            # component of bvec and setting it to 0.
            vec2 = list(bvec)
            minc = min(vec2)
            vec2[vec2.index(min(vec2))] = 0 if abs(minc) >= 1e-5 else 1
            # Get a vector orthogonal to bvec and vec2.
            a = normalize_vector(np.cross(bvec, vec2))
            return a

        else:
            return normalize_vector(self.centroid(conformer) -
                                    self.bonder_centroid(conformer))

    def fg_bonder_centroid(self, fg_id, conformer=-1):
        """
        The centroid of bonder atoms in a functional group.

        Parameters
        ----------
        fg_id : :class:`int`
            The id of the functional group.

        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        :class:`numpy.ndarray`
            The coordinates of a bonder centroid.

        """

        fg = self.func_groups[fg_id]
        bonder_coords = (self.atom_coords(bonder, conformer)
                         for bonder in fg.bonder_ids)
        return sum(bonder_coords) / len(fg.bonder_ids)

    def fg_distance(self, fg1, fg2, conformer=-1):
        """
        The distance between the bonder centroids of two fgs.

        Parameters
        ----------
        fg1 : :class:`int`
            The ``fg_id`` of the first fg.

        fg2 : :class:`int`
            The ``fg_id`` of the second fg.

        conformer : :class:`int`, optional
            The id of the conformer to use.

        Returns
        -------
        :class:`float`
            The distance between `fg1` and `fg2`.

        """

        c1 = self.fg_bonder_centroid(fg1, conformer)
        c2 = self.fg_bonder_centroid(fg2, conformer)
        return euclidean(c1, c2)

    def functional_groups(self, fg_names):
        """
        Finds given functional groups in the molecule.

        Parameters
        ----------
        fg_names : :class:`list` of :class:`str`
            The names of functional groups which are to be found
            within the molecule.

        Returns
        -------
        :class:`list` of :class:`.FunctionalGroup`
            A :class:`list` holding a :class:`.FunctionalGroup``
            instance for every matched functional group in the
            molecule.

        """

        # Ensure that given the same fg names in a different order,
        # the atoms get assigned to the a functional group with the
        # same id.
        fg_names = sorted(fg_names)

        func_groups = []
        for fg_name in fg_names:
            fg_info = fg_infos[fg_name]

            # Find all fg atoms.
            fg_query = rdkit.MolFromSmarts(fg_info.fg_smarts)
            fg_atoms = self.mol.GetSubstructMatches(fg_query)

            # Find all bonder atoms.
            bonder_atoms = [[] for i in range(len(fg_atoms))]

            for match in fg_info.bonder_smarts:
                query = rdkit.MolFromSmarts(match.smarts)
                atoms = set(flatten(
                    self.mol.GetSubstructMatches(query)
                ))

                # Get all the bonders grouped by fg.
                bonders = [[aid for aid in fg if aid in atoms]
                           for fg in fg_atoms]

                for fg_id, fg in enumerate(bonders):
                    bonder_atoms[fg_id].extend(fg[:match.n])

            # Find all deleter atoms.
            deleter_atoms = [[] for i in range(len(fg_atoms))]
            for match in fg_info.del_smarts:
                query = rdkit.MolFromSmarts(match.smarts)
                atoms = set(flatten(
                    self.mol.GetSubstructMatches(query)
                ))

                # Get all deleters grouped by fg.
                deleters = [[aid for aid in fg if aid in atoms]
                            for fg in fg_atoms]

                for fg_id, fg in enumerate(deleters):
                    deleter_atoms[fg_id].extend(fg[:match.n])

            for atom_ids in zip(fg_atoms, bonder_atoms, deleter_atoms):
                fg, bonders, deleters = atom_ids
                fg = FunctionalGroup(id_=len(func_groups),
                                     atom_ids=fg,
                                     bonder_ids=tuple(bonders),
                                     deleter_ids=tuple(deleters),
                                     info=fg_info)
                func_groups.append(fg)

        return func_groups

    def json(self, include_attrs=None):
        """
        Returns a JSON representation of the molecule.

        The representation has the following form:

        .. code-block:: python

            {
                'class' : 'StructUnit',
                'mol_block' : '''A string holding the V3000 mol
                                 block of the molecule.''',
                'note' : 'This molecule is nice.',
                'name' : 'benzene',
                'atom_props': {0: {'prop1': 1.0,
                                   'prop2': 'value1'}}
            }

        Parameters
        ----------
        include_attrs : :class:`list` of :class:`str`, optional
            The names of attributes of the molecule to be added to
            the JSON. Each attribute is saved as a string using
            :func:`repr`.

        Returns
        -------
        :class:`dict`
            A :class:`dict` which represents the molecule.

        """

        if include_attrs is None:
            include_attrs = []

        fg_names = [info.name for info in self.func_group_infos]
        conformers = [
            (conf.GetId(), self.mdl_mol_block(conformer=conf.GetId()))
            for conf in self.mol.GetConformers()
        ]
        json = {

            'class': self.__class__.__name__,
            'func_groups': fg_names,
            'conformers': conformers,
            'note': self.note,
            'name': self.name,
            'atom_props': self.atom_props

        }

        json.update(
            {attr: repr(getattr(self, attr)) for attr in include_attrs}
        )

        return json

    @classmethod
    def _json_init(cls, json_dict):
        """
        Completes a JSON initialization.

        This function is not to be used. Use :meth:`Molecule.load`
        for loading instances from a JSON string. That function will
        automatically call this one.

        Parameters
        ----------
        json_dict : :class:`dict`
            A dictionary holding the attribute data of the molecule.

        Returns
        -------
        None : :class:`NoneType`

        """

        d = dict(json_dict)
        d.pop('class')
        first_conf, *confs = d.pop('conformers')
        conf_id, mol_block = first_conf
        name = d.pop('name') if d.pop('load_names') else ''
        obj = cls(mol_block,
                  d.pop('func_groups'),
                  name,
                  d.pop('note'))

        obj.mol.GetConformer().SetId(conf_id)
        obj.atom_props = defaultdict(dict)
        obj.atom_props.update(d.pop('atom_props'))

        for conf_id, mol_block in confs:
            conf_mol = rdkit.MolFromMolBlock(molBlock=mol_block,
                                             removeHs=False,
                                             sanitize=False)
            conf = conf_mol.GetConformer()
            conf.SetId(conf_id)
            obj.mol.AddConformer(conf)

        for attr, val in d.items():
            setattr(obj, attr, eval(val))

        return obj

    @staticmethod
    def gen_key(rdkit_mol, functional_groups):
        """
        Generates the key used when caching the molecule.

        Parameters
        ----------
        rdkit_mol : :class:`rdkit.Mol`
            An ``rdkit`` instance of the molecule.

        functional_groups : :class:`tuple` of :class:`str`
            The name of the functional groups being used to make
            macromolecules.

        Returns
        -------
        :class:`tuple`
            The key used for caching the molecule. Has the form

            .. code-block:: python

                ('amine', 'bromine', 'InChIString')

        """

        functional_groups = sorted(functional_groups)
        return (*functional_groups, rdkit.MolToInchi(rdkit_mol))

    def minimize_theta(self, v1, v2, axis, centroid, conformer=-1):
        """
        Rotates the molecule to minimize angle between `v1` and `v2`.

        The rotation is done about the vector `axis`.

        Parameters
        ----------
        v1 : :class:`numpy.ndarray`
            The vector which is rotated.

        v2 : :class:`numpy.ndarray`
            The vector which is stationary.

        axis : :class:`numpy.ndarray`
            The vector about which the rotation happens.

        centroid : :class:`numpy.ndarray`
            The position vector at the center of the rotation.

        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        None : :class:`NoneType`

        """

        # If the vector being rotated is not finite exit. This is
        # probably due to a planar molecule.
        if not all(np.isfinite(x) for x in v1):
            return

        # Save the initial position and change the origin to
        # `centroid`.
        iposition = self.centroid(conformer)
        self.mol = self.shift(-centroid, conformer)

        # 1. First transform the problem.
        # 2. The rotation axis is set equal to the z-axis.
        # 3. Apply this transformation to all vectors in the problem.
        # 4. Take only the x and y components of `v1` and `v2`.
        # 5. Work out the angle between them.
        # 6. Apply that rotation along the original rotation axis.

        rotmat = rotation_matrix(axis, [0, 0, 1])
        tstart = np.dot(rotmat, v1)
        tstart = np.array([tstart[0], tstart[1], 0])

        # If the `tstart` vector is 0 after these transformations it
        # means that it is parallel to the rotation axis, stop.
        if np.allclose(tstart, [0, 0, 0], atol=1e-8):
            self.set_position(iposition, conformer)
            return

        tend = np.dot(rotmat, v2)
        tend = np.array([tend[0], tend[1], 0])
        angle = vector_theta(tstart, tend)

        # Check in which direction the rotation should go.
        # This is done by applying the rotation in each direction and
        # seeing which one leads to a smaller theta.
        r1 = rotation_matrix_arbitrary_axis(angle, [0, 0, 1])
        t1 = vector_theta(np.dot(r1, tstart), tend)
        r2 = rotation_matrix_arbitrary_axis(-angle, [0, 0, 1])
        t2 = vector_theta(np.dot(r2, tstart), tend)

        if t2 < t1:
            angle *= -1

        rot_mat = rotation_matrix_arbitrary_axis(angle, axis)
        pos_mat = self.mol.GetConformer(conformer).GetPositions().T
        new_pos_mat = np.dot(rot_mat, pos_mat)
        self.set_position_from_matrix(new_pos_mat, conformer)
        self.set_position(iposition, conformer)

    def rotate2(self, theta, axis, conformer=-1):
        """
        Rotates the molecule by `theta` about `axis`.

        The rotation occurs about the centroid of the bonder atoms.

        Parameters
        ----------
        theta : :class:`float`
            The size of the rotation in radians.

        axis : :class:`numpy.ndarray`
            The axis about which rotation happens.

        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        None : :class:`NoneType`

        """

        # Save the origin position of the bonder atom centroid.
        og_position = self.bonder_centroid(conformer)
        # Change the position of the centroid of the bonder atoms to
        # the origin so that the rotation occurs about this point.
        self.set_bonder_centroid([0, 0, 0], conformer)
        # Get the rotation matrix.
        rot_mat = rotation_matrix_arbitrary_axis(theta, axis)
        # Apply the rotation on the original atomic coordinates to get
        # the new ones.
        pos_mat = self.mol.GetConformer(conformer).GetPositions().T
        new_pos_mat = np.dot(rot_mat, pos_mat)
        # Set the atomic positions to the new coordinates.
        self.set_position_from_matrix(new_pos_mat, conformer)
        # Return the centroid to its original position.
        self.set_bonder_centroid(og_position, conformer)

    def set_bonder_centroid(self, position, conformer=-1):
        """
        Move the molecule so that the bonder centroid is on `position`.

        Parameters
        ----------
        position : :class:`numpy.ndarray`
            An array holding the desired the position. It holds
            the x, y and z coordinates, respectively.

        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        :class:`rdkit.Mol`
            The ``rdkit`` molecule after it has been shifted. The same
            instance as in :attr:`~Molecule.mol`.

        """

        conf_id = self.mol.GetConformer(conformer).GetId()

        center = self.bonder_centroid(conf_id)
        shift = position - center
        new_conf = self.shift(shift, conf_id).GetConformer()
        new_conf.SetId(conf_id)

        # Make sure the rkdit molecule has only one conformer.
        self.mol.RemoveConformer(conf_id)
        self.mol.AddConformer(new_conf)

        return self.mol

    def set_orientation2(self, end, conformer=-1):
        """
        Rotate the molecule so that bonder atoms lie on `end`.

        The molecule is rotated about the centroid of the bonder atoms.
        It is rotated so that the direction vector running between the
        2 bonder centroids is aligned with the vector `end`.

        Parameters
        ----------
        end : :class:`numpy.ndarray`
            The vector with which the molecule's bonder atoms should be
            aligned.

        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        :class:`rdkit.Mol`
            The ``rdkit`` molecule in :attr:`~Molecule.mol`.

        """

        start = self.centroid() - self.bonder_centroid()
        return self._set_orientation2(start, end, conformer)

    def _set_orientation2(self, start, end, conformer):
        """
        Rotates the molecule by a rotation from `start` to `end`.

        Given two direction vectors, `start` and `end`, this method
        applies the rotation required transform `start` to `end` on
        the molecule. The rotation occurs about the centroid of the
        molecule.

        For example, if the `start` and `end` vectors
        are 45 degrees apart, a 45 degree rotation will be applied to
        the molecule. The rotation will be along the appropriate
        direction.

        The great thing about this method is that you as long as you
        can associate a gemotric feature of the molecule with a vector,
        then the molecule can be roatated so that this vector is
        aligned with `end`. The defined vector can be virtually
        anything. This means that any geometric feature of the molecule
        can be easily aligned with any arbitrary axis.

        Notes
        -----
        The difference between this method and
        :meth:`~Molecule.set_orientation` is about which point the
        rotation occurs: centroid of bonder atoms versus centroid of
        entire molecule, respectively.

        Parameters
        ----------
        start : :class:`numpy.ndarray`
            A vector which is to be rotated so that it transforms to
            the `end` vector.

        end : :class:`numpy.ndarray`
            This array holds the vector, onto which `start` is rotated.

        conformer : :class:`int`
            The conformer to use.

        Returns
        -------
        :class:`rdkit.Mol`
            The ``rdkit`` molecule in :attr:`~Molecule.mol`.

        """

        # Normalize the input direction vectors.
        start = normalize_vector(start)
        end = normalize_vector(end)

        # Record the position of the molecule then translate the bonder
        # atom centroid to the origin. This is so that the rotation
        # occurs about this point.
        og_center = self.bonder_centroid(conformer)
        self.set_bonder_centroid(np.array([0, 0, 0]), conformer)

        # Get the rotation matrix.
        rot_mat = rotation_matrix(start, end)

        # Apply the rotation matrix to the atomic positions to yield
        # the new atomic positions.
        pos_mat = self.mol.GetConformer(conformer).GetPositions().T
        new_pos_mat = np.dot(rot_mat, pos_mat)

        # Set the positions in the rdkit molecule.
        self.set_position_from_matrix(new_pos_mat, conformer)
        self.set_bonder_centroid(og_center, conformer)

        return self.mol

    def shift_fgs(self, new_ids, shift):
        """
        Yield new functional groups with atomic ids shifted.

        Parameters
        ----------
        new_ids : :class:`iterable` of :class:`int`
            The ids assigned to the new functional groups.

        num_atoms : :class:`int`
            The number to shift each atom id by.

        Yields
        ------
        :class:`.FunctionalGroup`
            A functional group from :attr:`~StructUnit.func_groups`
            with atomic ids shifted by `shift`.

        """

        for id_, fg in zip(new_ids, self.func_groups):
            yield fg.shifted_fg(id_=id_, shift=shift)

    def similar_molecules(self, mols):
        """
        Returns molecules from `mols` ordered by similarity.

        The most similar molecule is at index 0.

        This method uses the Morgan fingerprints of radius 4 to
        evaluate how similar the molecules in `mols` are.

        Parameters
        ----------
        mols : :class:`iterable` of :class:`rdkit.Mol`
            A group of molecules to which similarity is compared.

        Returns
        -------
        :class:`list`
            A :class:`list` of the form,

            .. code-block:: python

                returned_list = [(8.9, mol1), (7.3, mol2), (3.4, mol3)]

            where the :class:`float` is the similarity of a given
            molecule in `mols` while the ```mol`` is corresponding
            ``rdkit`` molecule. Most similar molecule yielded first.

        """

        # First get the fingerprint of `self`.
        rdkit.GetSSSR(self.mol)
        self.mol.UpdatePropertyCache(strict=False)
        fp = rdkit.GetMorganFingerprint(self.mol, 4)

        # For every structure file in the database create a rdkit
        # molecule. Place these in a list.
        similarities = []
        for mol in mols:
            rdkit.GetSSSR(mol)
            mol.UpdatePropertyCache(strict=False)
            mol_fp = rdkit.GetMorganFingerprint(mol, 4)
            similarity = DataStructs.DiceSimilarity(fp, mol_fp)
            similarities.append((similarity, mol))

        return sorted(similarities, reverse=True, key=lambda x: x[0])

    @classmethod
    def smiles_init(cls,
                    smiles,
                    functional_groups=None,
                    random_seed=4,
                    note="",
                    name=""):
        """
        Initialize from a SMILES string.

        The structure of the molecule is embedded using
        :func:`rdkit.ETKDG()`.

        Parameters
        ----------
        smiles : :class:`str`
            A SMILES string of the molecule.

        functional_groups : :class:`list` of :class:`str`, optional
            The name of the functional groups which are to have atoms
            tagged. If no functional groups are provided to this
            parameter, no tagging is done.

        random_seed : :class:`int`, optional
            Random seed passed to :func:`rdkit.ETKDG`

        note : :class:`str`, optional
            A note or comment about the molecule.

        name : :class:`str`, optional
            A name which can be optionally given to the molcule for
            easy identification.

        Returns
        -------
        :class:`StructUnit`
            A :class:`StructUnit` instance of the molecule in `smarts`.

        """

        if functional_groups is None:
            functional_groups = ()

        mol = rdkit.MolFromSmiles(smiles)
        rdkit.SanitizeMol(mol)
        H_mol = rdkit.AddHs(mol)
        key = cls.gen_key(H_mol, functional_groups)
        if key in cls.cache and OPTIONS['cache']:
            return cls.cache[key]

        params = rdkit.ETKDGv2()
        params.randomSeed = random_seed
        for i in range(100):
            failed = rdkit.EmbedMolecule(mol, params) == -1
            if failed:
                params.randomSeed += 1
            else:
                break

        if params.randomSeed != random_seed:
            msg = ('Embedding with seed value of '
                   f'"{random_seed}" failed. Using alternative value'
                   f' of "{params.randomSeed}" was successful.')
            logger.warning(msg)

        mol = rdkit.AddHs(mol, addCoords=True)
        mol.GetConformer()
        obj = cls.__new__(cls)
        obj.file = smiles
        obj.key = key
        obj.mol = mol
        obj.func_groups = tuple(
            obj.functional_groups(functional_groups)
        )
        obj.func_group_infos = tuple(dedupe(
            iterable=(fg.info for fg in obj.func_groups),
            key=lambda info: info.name
        ))

        Molecule.__init__(obj, name, note)

        cls.cache[key] = obj
        return obj

    def __str__(self):
        return "{} {}".format(self.__class__.__name__, list(self.key))

    def __repr__(self):
        return str(self)


class StructUnit2(StructUnit):
    """
    Represents building blocks with 2 functional groups.

    """

    def set_orientation2(self, end, conformer=-1):
        """
        Rotate the molecule so that bonder atoms lie on `end`.

        The molecule is rotated about the centroid of the bonder atoms.
        It is rotated so that the direction vector running between the
        2 bonder centroids is aligned with the vector `end`.

        Parameters
        ----------
        end : :class:`numpy.ndarray`
            The vector with which the molecule's bonder atoms should be
            aligned.

        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        :class:`rdkit.Mol`
            The ``rdkit`` molecule in :attr:`~Molecule.mol`.

        """

        *_, start = next(self.bonder_direction_vectors(conformer))
        return self._set_orientation2(start, end, conformer)

    def minimize_theta2(self, vector, axis, conformer=-1):
        """
        Rotates molecule about `axis` to minimze theta with `vector`.

        The molecule is rotated about `axis` passing through the bonder
        centroid. It is rotated so that the vector between the bonder
        and molecular centroids lies on the same plane as `vector`.

        Parameters
        ----------
        vector : :class:`numpy.ndarray`
            The vector to which the distance should be minimized.

        axis : :class:`numpy.ndarray`
            The direction vector along which the rotation happens.

        conformer : :class:`int`, optional
            The id of the conformer to be used.

        Returns
        -------
        None : :class:`NoneType`

        """

        self.minimize_theta(
            v1=self.centroid_centroid_dir_vector(conformer),
            v2=vector,
            axis=axis,
            centroid=self.bonder_centroid(conformer),
            conformer=conformer)


class StructUnit3(StructUnit):
    """
    Represents building blocks with 3 or more functional groups.

    """

    def bonder_plane(self, conformer=-1):
        """
        Returns coeffs of the plane formed by the bonder centroids.

        A plane is defined by the scalar plane equation::

            ax + by + cz = d.

        This method returns the ``a``, ``b``, ``c`` and ``d``
        coefficients of this equation for the plane formed by the
        bonder centroids. The coefficents ``a``, ``b`` and ``c``
        describe the normal vector to the plane. The coefficent ``d``
        is found by substituting these coefficients along with the
        x, y and z variables in the scalar equation and
        solving for ``d``. The variables x, y and z are
        substituted by the coordinates of some point on the plane. For
        example, the position of one of the bonder centroids.

        Parameters
        ----------
        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        :class:`numpy.ndarray`
            This array has the form ``[a, b, c, d]`` and represents the
            scalar equation of the plane formed by the bonder
            centroids.

        References
        ----------
        https://tinyurl.com/okpqv6

        """

        centroid = next(self.bonder_centroids(conformer))
        d = -np.sum(self.bonder_plane_normal(conformer) * centroid)
        return np.append(self.bonder_plane_normal(conformer), d)

    def bonder_plane_normal(self, conformer=-1):
        """
        Returns the normal to the plane formed by bonder centroids.

        The normal of the plane is defined such that it goes in the
        direction toward the centroid of the molecule.

        Parameters
        ----------
        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        :class:`numpy.ndarray`
            A unit vector which describes the normal to the plane of
            the bonder centroids.

        """

        if len(self.func_groups) < 3:
            raise ValueError(("StructUnit3 molecule "
                             "has fewer than 3 functional groups."))

        direction_vectors = (
            v for *_, v in self.bonder_direction_vectors(conformer)
        )
        v1, v2 = it.islice(direction_vectors, 2)

        normal_v = normalize_vector(np.cross(v1, v2))

        theta = vector_theta(
                    normal_v,
                    self.centroid_centroid_dir_vector(conformer))

        if theta > np.pi/2:
            normal_v *= -1

        return normal_v

    def minimize_theta2(self, fg_id, vector, axis, conformer=-1):
        """
        Rotate molecule to minimize angle between `fg_id` and `vector`.

        The rotation is done about `axis` and is centered on the
        bonder centroid, as given by
        :meth:`~.StructUnit.bonder_centroid`.

        Parameters
        ----------
        fg_id : :class:`int`
            The id of functional group which is to have angle
            minimized.

        vector : :class:`numpy.ndarray`
            A vector with which the angle is minimized.

        axis : :class:`numpy.ndarray`
            The vector about which the rotation happens.

        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        None : :class:`NoneType`

        """

        centroid = self.bonder_centroid(conformer)
        fg = self.func_groups[fg_id]

        bonder_centroid = self.atom_centroid(fg.bonder_ids)
        v1 = bonder_centroid - centroid
        self.minimize_theta(v1, vector, axis, centroid, conformer)

    def set_orientation2(self, end, conformer=-1):
        """
        Rotates the molecule so the plane normal is aligned with `end`.

        Here "plane normal" referes to the normal of the plane formed
        by the bonder centroids. The molecule is rotated about
        :meth:`~Molecule.bonder_centroid`. The rotation results in the
        normal of the plane being aligned with `end`.

        Parameters
        ----------
        end : :class:`numpy.ndarray`
            The vector with which the normal of plane of bonder
            centroids shoould be aligned.

        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        :class:`rdkit.Mol`
            The ``rdkit`` molecule in :attr:`~Molecule.mol`.

        """

        start = self.bonder_plane_normal(conformer)
        return self._set_orientation2(start, end, conformer)
