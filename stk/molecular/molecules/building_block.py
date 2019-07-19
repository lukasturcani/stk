"""
Defines the :class:`BuildingBlock` class.

"""


import logging
import os
import numpy as np
import itertools as it
import rdkit.Chem.AllChem as rdkit
from glob import glob
from functools import partial
from scipy.spatial.distance import euclidean

from . import elements
from .elements import Atom, Bond
from .molecule import Molecule
from ..functional_groups import FunctionalGroup
from ..functional_groups import functional_group_infos as fg_infos
from ..functional_groups import functional_groups as fgs
from ...utilities import (
    flatten,
    normalize_vector,
    vector_theta,
    remake,
    dedupe
)


logger = logging.getLogger(__name__)


class BuildingBlock(Molecule):
    """
    Represents the building blocks of a :class:`.ConstructedMolecule`.

    :class:`BuildingBlock` represents building block molecules, which
    are either entire molecules or molecular fragments used for the
    construction of a :class:`.ConstructedMolecule`.

    Attributes
    ----------
    _init_funcs : :class:`dict`
        This is a class attribute. It maps file extensions to
        functions which can be used to initialize
        :class:`BuildingBlock` from that file type.

    func_groups : :class:`tuple` of :class:`.FunctionalGroup`
        The functional groups present in the molecule. The
        :attr:`.FunctionalGroup.id` is equal to the index.

    _key : :class:`object`
        Extends :class:`.Molecule._key`. :class:`BuildingBlock`
        molecules with the same structure and with the same functional
        groups will have equal keys. Used for caching.

    Methods
    -------
    :meth:`__init__`
    :meth:`init_from_random_file`
    :meth:`init_from_smiles`
    :meth:`get_bonder_ids`
    :meth:`get_bonder_centroids`
    :meth:`get_bonder_plane`
    :meth:`get_bonder_plane_normal`
    :meth:`get_bonder_distances`
    :meth:`get_bonder_direction_vectors`
    :meth:`get_centroid_centroid_direction_vector`
    :meth:`get_functional_groups`
    :meth:`to_dict`
    :meth:`shift_fgs`

    """

    _init_funcs = {
        '.mol': partial(
            rdkit.MolFromMolFile,
            sanitize=False,
            removeHs=False
        ),

        '.sdf': partial(
            rdkit.MolFromMolFile,
            sanitize=False,
            removeHs=False
        ),

        '.mol2': partial(
            rdkit.MolFromMol2File,
            sanitize=False,
            removeHs=False
        ),

        '.pdb': partial(
            rdkit.MolFromPDBFile,
            sanitize=False,
            removeHs=False
        ),
    }

    def __init__(self, mol, functional_groups=None, use_cache=False):
        """
        Initialize a :class:`BuildingBlock` instance.

        Parameters
        ----------
        mol : :class:`str` or :class:`rdkit.Mol`
            Can be one of 2 things:

                1. A path to a molecular structure file.
                2. A :class:`rdkit.Mol` object.

            Supported file types are:

                #. ``.mol``, ``.sdf`` - MDL V3000 MOL file
                #. ``.pdb`` - PDB file
                #. ``.mol2`` - SYBYL MOL2 file

        functional_groups : :class:`list` of :class:`str`, optional
            The names of the functional groups which are to have atoms
            tagged. If ``None``, a functional group name found in the
            path `file`  is used. If no functional groups are provided
            to this parameter and the name of one is not present in
            `file`, no tagging is done.

        use_cache : :class:`bool`, optional
            If ``True``, a new :class:`.BuildingBlock` will
            not be made if a cached and identical one already exists,
            the one which already exists will be returned. If ``True``
            and a cached, identical :class:`BuildingBlock` does not
            yet exist the created one will be added to the cache.

        Raises
        ------
        :class:`TypeError`
            If the file type cannot be used for initialization.

        """

        if isinstance(mol, str):
            if os.path.exists(mol):
                _, ext = os.path.splitext(mol)

                if ext not in self._init_funcs:
                    raise TypeError(
                        f'Unable to initialize from "{ext}" files.'
                    )
                rdkit_mol = remake(self._init_funcs[ext](mol))

            else:
                rdkit_mol = remake(
                    rdkit.MolFromMolBlock(
                        molBlock=mol,
                        removeHs=False,
                        sanitize=False
                    )
                )

        elif isinstance(mol, rdkit.Mol):
            rdkit_mol = remake(mol)

        # Update the property cache of each atom. This updates things
        # like valence.
        for atom in rdkit_mol.GetAtoms():
            atom.UpdatePropertyCache()

        atoms = tuple(
            Atom(a.GetIdx(), a.GetAtomicNum())
            for a in rdkit_mol.GetAtoms()
        )
        bonds = tuple(
            Bond(
                atoms[b.GetBeginAtomIdx()],
                atoms[b.GetEndAtomIdx()],
                b.GetBondTypeAsDouble()
            )
            for b in rdkit_mol.GetBonds()
        )

        position_matrix = rdkit_mol.GetConformer().GetPositions()
        super().__init__(atoms, bonds, position_matrix)

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
            self.get_functional_groups(functional_groups)
        )

    @classmethod
    def init_from_random_file(
        cls,
        file_glob,
        functional_groups=None,
        use_cache=False
    ):
        """
        Pick a random file from `file_glob` to initialize from.

        Parameters
        ----------
        file_glob : :class:`str`
            A glob specifying files, one of which is used to initialize
            a :class:`.BuildingBlock` at random.

        functional_groups : :class`list` of :class:`str`, optional
            The names of a functional groups which the molecule
            selected from the `file_glob` should have.

        use_cache : :class:`bool`, optional
            If ``True``, a new :class:`.BuildingBlock` will
            not be made if a cached and identical one already exists,
            the one which already exists will be returned. If ``True``
            and a cached, identical :class:`BuildingBlock` does not
            yet exist the created one will be added to the cache.

        Returns
        -------
        :class:`BuildingBlock`
            A random molecule from `file_glob`.

        Raises
        ------
        :class:`RuntimeError`
            If no files in `file_glob` could be initialized from.

        """

        files = glob(file_glob)
        np.random.shuffle(file_glob)

        for mol_file in files:
            try:
                return cls(mol_file, functional_groups, use_cache)

            except Exception:
                msg = (
                    'Could not initialize '
                    f'{cls.__name__} from {mol_file}.'
                )
                logger.warning(msg)
        raise RuntimeError(
            f'No files in "{file_glob}" could be initialized from.'
        )

    @classmethod
    def init_from_smiles(
        cls,
        smiles,
        functional_groups=None,
        random_seed=4,
        use_cache=False
    ):
        """
        Initialize from a SMILES string.

        The structure of the molecule is embedded using
        :func:`rdkit.ETKDGv2()`.

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

        use_cache : :class:`bool`, optional
            If ``True``, a new :class:`.BuildingBlock` will
            not be made if a cached and identical one already exists,
            the one which already exists will be returned. If ``True``
            and a cached, identical :class:`BuildingBlock` does not
            yet exist the created one will be added to the cache.

        Returns
        -------
        :class:`BuildingBlock`
            A :class:`BuildingBlock` instance of the molecule in
            `smiles`.

        """

        if functional_groups is None:
            functional_groups = ()

        mol = rdkit.AddHs(rdkit.MolFromSmiles(smiles))

        params = rdkit.ETKDGv2()
        params.randomSeed = random_seed
        for i in range(100):
            failed = rdkit.EmbedMolecule(mol, params) == -1
            if failed:
                params.randomSeed += 1
            else:
                break

        if params.randomSeed != random_seed:
            msg = (
                'Embedding with seed value of '
                f'"{random_seed}" failed. Using alternative value '
                f'of "{params.randomSeed}" was successful.'
            )
            logger.warning(msg)

        return cls(
            mol=mol,
            functional_groups=functional_groups,
            use_cache=use_cache
        )

    @classmethod
    def _init_from_dict(cls, mol_dict, use_cache):
        """
        Intialize from a :class:`dict` representation.

        Parameters
        ----------
        mol_dict : :class:`dict`
            A :class:`dict` representation of a molecule generated
            by :meth:`to_dict`.

        use_cache : :class:`bool`
            If ``True``, a new instance will not be made if a cached
            and identical one already exists, the one which already
            exists will be returned. If ``True`` and a cached,
            identical instance does not yet exist the created one will
            be added to the cache.

        Returns
        -------
        :class:`BuildingBlock`
            The molecule described by `mol_dict`.

        """

        d = dict(mol_dict)
        d.pop('class')
        mol_block = d.pop('mol_block')
        rdkit_mol = remake(rdkit.MolFromMolBlock(
            molBlock=mol_block,
            removeHs=False,
            sanitize=False
        ))
        obj = cls(rdkit_mol, d.pop('func_groups'), use_cache)
        obj.atoms = eval(d.pop('atoms'), vars(elements))
        for attr, val in d.items():
            setattr(obj, attr, eval(val))
        return obj

    def get_bonder_ids(self, fg_ids=None):
        """
        Yield ids of bonder atoms.

        Parameters
        ----------
        fg_ids : :class:`iterable` of :class:`int`
            The ids of functional groups whose bonder atoms should be
            yielded. If ``None`` then all bonder atom ids in the
            :class:`.BuildingBlock` will be yielded.

        Yields
        ------
        :class:`int`
            The id of a bonder atom.

        """

        if fg_ids is None:
            fg_ids = range(len(self.func_groups))

        for fg_id in fg_ids:
            for atom_id in self.func_groups[fg_id].bonder_ids:
                yield atom_id

    def get_bonder_centroids(self, fg_ids=None):
        """
        Yield the centroids of bonder atoms in functional groups.

        A bonder centroid is the centroid of all bonder atoms in a
        particular functional group.

        The bonder centroids are yielded in the order given by
        `fg_ids`.

        Parameters
        ----------
        fg_ids : :class:`iterable` of :class:`int`
            The ids of functional groups to be used.
            If ``None`` then all functional groups are used and
            centroids are yielded in ascending order of functional
            group id.

        Yields
        ------
        :class:`numpy.ndarray`
            The centroid of a functional groups

        """

        if fg_ids is None:
            fg_ids = range(len(self.func_groups))

        for fg_id in fg_ids:
            yield self.get_centroid(
                atom_ids=self.func_groups[fg_id].bonder_ids
            )

    def get_bonder_plane(self, fg_ids=None):
        """
        Return coeffs of the plane formed by the bonder centroids.

        A bonder centroid is the centroid of all bonder atoms in a
        particular functional group.

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
        fg_ids : :class:`iterable` of :class:`int`
            The ids of functional groups used to contruct the plane.
            If there are more than three, a plane of best fit through
            the bonder centroids of the functional groups will be made.
            If ``None``, all functional groups in the
            :class:`BuildingBlock` will be used.

        Returns
        -------
        :class:`numpy.ndarray`
            This array has the form ``[a, b, c, d]`` and represents the
            scalar equation of the plane formed by the bonder
            centroids.

        Raises
        ------
        :class:`RuntimeError`
            If there are not at least 3 functional groups, which is
            necessary to define a plane.

        References
        ----------
        https://tinyurl.com/okpqv6

        """

        if fg_ids is None:
            fg_ids = range(len(self.func_groups))
        else:
            # Iterable is used multiple times.
            fg_ids = list(fg_ids)

        centroid = self.get_centroid(
            atom_ids=self.func_groups[fg_ids[0]].bonder_ids
        )
        normal = self.get_bonder_plane_normal(
            fg_ids=fg_ids
        )
        d = np.sum(normal * centroid)
        return np.append(normal, d)

    def get_bonder_plane_normal(self, fg_ids=None):
        """
        Return the normal to the plane formed by bonder centroids.

        A bonder centroid is the centroid of all bonder atoms in a
        particular functional group.

        The normal of the plane is defined such that it goes in the
        direction toward the centroid of the molecule.

        Parameters
        ----------
        fg_ids  : :class:`iterable` of :class:`int`, optional
            The ids of functional groups used to contruct the plane.
            If there are more than three, a plane of best fit through
            the bonder centroids of the functional groups will be made.
            If ``None``, all functional groups in the
            :class:`BuildingBlock` will be used.

        Returns
        -------
        :class:`numpy.ndarray`
            A unit vector which describes the normal to the plane of
            the bonder centroids.

        Raises
        ------
        :class:`RuntimeError`
            If there are not at least 3 functional groups, which is
            necessary to define a plane.

        """

        if fg_ids is None:
            fg_ids = range(len(self.func_groups))
        else:
            # The iterable is used mutliple times.
            fg_ids = list(fg_ids)

        if len(fg_ids) < 3:
            raise RuntimeError(
                'At least 3 functional groups '
                'are necessary to create a plane.'
            )

        centroids = np.array(list(self.get_bonder_centroids(
            fg_ids=fg_ids
        )))
        bonder_centroid = self.get_centroid(
            atom_ids=self.get_bonder_ids(fg_ids=fg_ids)
        )
        normal = np.linalg.svd(centroids - bonder_centroid)[-1][2, :]
        cc_vector = self.get_centroid_centroid_direction_vector(
            fg_ids=fg_ids
        )
        if vector_theta(normal, cc_vector) > np.pi/2:
            normal *= -1
        return normal

    def get_bonder_distances(self, fg_ids=None):
        """
        Yield distances between pairs of bonder centroids.

        A bonder centroid is the centroid of all bonder atoms in a
        particular functional group.

        Parameters
        ----------
        fg_ids : :class:`iterable` of :class:`int`
            The ids of functional groups to be used.
            If ``None`` then all functional groups are used.

        Yields
        ------
        :class:`tuple`
            A :class:`tuple` of the form ``(3, 54, 12.54)``.
            The first two elements are the ids of the involved
            functional groups and the third element is the distance
            between them.

        """

        if fg_ids is None:
            fg_ids = range(len(self.func_groups))

        # Iterator yielding tuples of form (fg_id, bonder_centroid)
        centroids = ((
            i, self.get_centroid(
                atom_ids=self.func_groups[i].bonder_ids
            ))
            for i in fg_ids
        )
        pairs = it.combinations(iterable=centroids, r=2)
        for (id1, c1), (id2, c2) in pairs:
            yield id1, id2, float(euclidean(c1, c2))

    def get_bonder_direction_vectors(
        self,
        fg_ids=None
    ):
        """
        Yield the direction vectors between bonder centroids.

        A bonder centroid is the centroid of all bonder atoms in a
        particular functional group.

        Parameters
        ----------
        fg_ids : :class:`iterable` of :class:`int`
            The ids of functional groups to be used.
            If ``None`` then all functional groups are used.

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

        if fg_ids is None:
            fg_ids = range(len(self.func_groups))

        # Iterator yielding tuples of form (fg_id, bonder_centroid)
        centroids = ((
            i, self.get_centroid(
                atom_ids=self.func_groups[i].bonder_ids
            ))
            for i in fg_ids
        )
        pairs = it.combinations(iterable=centroids, r=2)
        for (id1, c1), (id2, c2) in pairs:
            yield id2, id1, normalize_vector(c1-c2)

    def get_centroid_centroid_direction_vector(
        self,
        fg_ids=None
    ):
        """
        Return the direction vector between the 2 molecular centroids.

        The first molecular centroid is the centroid of the entire
        molecule. The second molecular centroid is the of the
        bonder atoms in the molecule.

        Parameters
        ----------
        fg_ids : :class:`iterable` of :class:`int`
            The ids of functional groups to be used for calculating the
            bonder centroid. If ``None`` then all functional groups are
            used.

        Returns
        -------
        :class:`numpy.ndarray`
            The normalized direction vector running from the centroid
            of the bonder atoms to the molecular centroid.

        """

        if fg_ids is None:
            fg_ids = range(len(self.func_groups))
        else:
            # This iterable gets used more than once.
            fg_ids = list(fg_ids)

        bonder_centroid = self.get_centroid(
            atom_ids=self.get_bonder_ids(fg_ids=fg_ids)
        )
        centroid = self.get_centroid()
        # If the bonder centroid and centroid are in the same position,
        # the centroid - centroid vector should be orthogonal to the
        # bonder direction vector.
        if np.allclose(centroid, bonder_centroid, 1e-5):
            *_, bvec = self.get_bonder_direction_vectors(
                fg_ids=fg_ids
            )
            # Construct a secondary vector by finding the minimum
            # component of bvec and setting it to 0.
            vec2 = list(bvec)
            minc = min(vec2)
            vec2[vec2.index(min(vec2))] = 0 if abs(minc) >= 1e-5 else 1
            # Get a vector orthogonal to bvec and vec2.
            return normalize_vector(np.cross(bvec, vec2))

        else:
            return normalize_vector(centroid - bonder_centroid)

    def get_functional_groups(self, fg_names):
        """
        Find functional groups in the molecule.

        Parameters
        ----------
        fg_names : :class:`list` of :class:`str`
            The names of the functional groups which are to be found
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

        mol = self.to_rdkit_mol()
        rdkit.SanitizeMol(mol)
        func_groups = []
        for fg_name in fg_names:
            fg_info = fg_infos[fg_name]

            # Find all fg atoms.
            fg_query = rdkit.MolFromSmarts(fg_info.fg_smarts)
            fg_atoms = mol.GetSubstructMatches(fg_query)

            # Find all bonder atoms.
            bonder_atoms = [[] for i in range(len(fg_atoms))]

            for match in fg_info.bonder_smarts:
                query = rdkit.MolFromSmarts(match.smarts)
                atoms = set(flatten(mol.GetSubstructMatches(query)))

                # Get all the bonders grouped by fg.
                bonders = [
                    [aid for aid in fg if aid in atoms]
                    for fg in fg_atoms
                ]

                for fg_id, fg in enumerate(bonders):
                    bonder_atoms[fg_id].extend(fg[:match.n])

            # Find all deleter atoms.
            deleter_atoms = [[] for i in range(len(fg_atoms))]
            for match in fg_info.del_smarts:
                query = rdkit.MolFromSmarts(match.smarts)
                atoms = set(flatten(mol.GetSubstructMatches(query)))

                # Get all deleters grouped by fg.
                deleters = [
                    [aid for aid in fg if aid in atoms]
                    for fg in fg_atoms
                ]

                for fg_id, fg in enumerate(deleters):
                    deleter_atoms[fg_id].extend(fg[:match.n])

            for atom_ids in zip(fg_atoms, bonder_atoms, deleter_atoms):
                fg, bonders, deleters = atom_ids
                fg = FunctionalGroup(
                    id=len(func_groups),
                    atom_ids=fg,
                    bonder_ids=tuple(bonders),
                    deleter_ids=tuple(deleters),
                    info=fg_info
                )
                func_groups.append(fg)

        return func_groups

    def to_dict(self, include_attrs=None):
        """
        Return a :class:`dict` representation of the molecule.

        Parameters
        ----------
        include_attrs : :class:`list` of :class:`str`, optional
            The names of additional attributes of the molecule to be
            added to the :class:`dict`. Each attribute is saved as a
            string using :func:`repr`.

        Returns
        -------
        :class:`dict`
            A :class:`dict` which represents the molecule. It has
            the form

            .. code-block:: python

                {
                    'class' : 'BuildingBlock',
                    'func_groups': ['amine', 'aldehyde'],
                    'mol_block' : '''A string holding the V3000 mol
                                     block of the molecule.''',
                    'atoms': [H(0), N(1), ...],
                }

        """

        if include_attrs is None:
            include_attrs = []

        fgs = list(dedupe(fg.info.name for fg in self.func_groups))
        d = {
            'class': self.__class__.__name__,
            'func_groups': fgs,
            'mol_block': self._to_mdl_mol_block(),
            'atoms': repr(self.atoms)
        }

        d.update(
            {attr: repr(getattr(self, attr)) for attr in include_attrs}
        )

        return d

    @staticmethod
    def _get_key(self, mol, functional_groups, use_cache):
        """
        Get the key used for caching the molecule.

        Parameters
        ----------
        mol : :class:`str` or :class:`rdkit.Mol`
            Can be one of 3 things:

                1. A path to a molecular structure file.
                2. A :class:`rdkit.Mol` object.
                3. V3000 MDL Mol block.

        functional_groups : :class:`list` of :class:`str`
            The names of the functional groups which are to have atoms
            tagged. If ``None``, a functional group name found in the
            path `file`  is used. If no functional groups are provided
            to this parameter and the name of one is not present in
            `file`, no tagging is done.

        use_cache : :class:`bool`
            This argument is ignored but included to be maintain
            compatiblity the the :meth:`__init__` signature.

        Returns
        -------
        :class:`tuple`
            The key used for caching the molecule. Has the form

            .. code-block:: python

                ('amine', 'bromine', 'InChIString')

        Raises
        ------
        :class:`TypeError`
            If the there is file type cannot be used for
            initialization.

        """

        if isinstance(mol, str):
            if os.path.exists(mol):
                _, ext = os.path.splitext(mol)

                if ext not in self._init_funcs:
                    raise TypeError(
                        f'Unable to initialize from "{ext}" files.'
                    )
                mol = remake(self._init_funcs[ext](mol))

            else:
                mol = remake(
                    rdkit.MolFromMolBlock(
                        molBlock=mol,
                        removeHs=False,
                        sanitize=False
                    )
                )

        elif isinstance(mol, rdkit.Mol):
            mol = remake(mol)

        if functional_groups is None:
            if isinstance(mol, str) and os.path.exists(mol):
                functional_groups = [
                    fg.name for fg in fgs if fg.name in mol
                ]
            else:
                functional_groups = []

        functional_groups = sorted(functional_groups)
        return (*functional_groups, rdkit.MolToInchi(mol))

    def shift_fgs(self, new_ids, shift):
        """
        Yield new functional groups with atomic ids shifted.

        Parameters
        ----------
        new_ids : :class:`iterable` of :class:`int`
            The ids assigned to the new functional groups.

        shift : :class:`int`
            The number to shift each atom id by.

        Yields
        ------
        :class:`.FunctionalGroup`
            A functional group from :attr:`~BuildingBlock.func_groups`
            with atomic ids shifted by `shift`.

        """

        for id, fg in zip(new_ids, self.func_groups):
            yield fg.shifted_fg(id=id, shift=shift)

    def __str__(self):
        return f'{self.__class__.__name__} {list(self._key)}'

    def __repr__(self):
        return str(self)
