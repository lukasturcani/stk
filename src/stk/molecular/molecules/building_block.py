"""
Building Block
==============

"""


import logging
import os
import numpy as np
import itertools as it
import rdkit.Chem.AllChem as rdkit
from glob import glob
from functools import partial
from scipy.spatial.distance import euclidean

from .. import atoms
from ..atoms import Atom
from ..bond import Bond
from .molecule import Molecule, _Molecule
from ...utilities import vector_angle, dedupe, remake


logger = logging.getLogger(__name__)


class BuildingBlock(Molecule, _Molecule):
    """
    Represents a building block of a :class:`.ConstructedMolecule`.

    A :class:`BuildingBlock` can represent either an entire molecule or
    a molecular fragments used to construct a
    :class:`.ConstructedMolecule`. The building block uses
    :class:`.FunctionalGroup` instances to identify which atoms are
    modified during construction.

    """

    # Maps file extensions to functions which can be used to
    # create an rdkit molecule from that file type.
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

        '.pdb': partial(
            rdkit.MolFromPDBFile,
            sanitize=False,
            removeHs=False
        ),
    }

    def __init__(self, smiles, functional_groups=None, random_seed=4):
        """
        Initialize a :class:`.BuildingBlock`.

        Notes
        -----
        The molecule is given 3D coordinates with
        :func:`rdkit.ETKDGv2()`.

        Parameters
        ----------
        smiles : :class:`str`
            A SMILES string of the molecule.

        functional_groups : :class:`iterable`, optional
            An :class:`iterable` of :class:`.FunctionalGroup` or
            :class:`.FunctionalGroupFactory` or both.
            :class:`.FunctionalGroup` instances are added to the
            building block and :class:`.FunctionalGroupFactory`
            instances are used to create :class:`.FunctionalGroup`
            instances the building block should hold.
            :class:`.FunctionalGroup` instances are used to identify
            which atoms are modified during
            :class:`.ConstructedMolecule` construction.

        random_seed : :class:`int`, optional
            Random seed passed to :func:`rdkit.ETKDGv2`

        Raises
        ------
        :class:`RuntimeError`
            If embedding the molecule fails.

        """

        if functional_groups is None:
            functional_groups = ()

        molecule = rdkit.AddHs(rdkit.MolFromSmiles(smiles))
        params = rdkit.ETKDGv2()
        params.randomSeed = random_seed
        if rdkit.EmbedMolecule(molecule, params) == -1:
            raise RuntimeError(
                f'Embedding with seed value of {random_seed} failed.'
            )
        identity_key = self._get_identity_key_from_rdkit_mol(
            molecule=molecule,
            functional_groups=functional_groups,
        )

        rdkit.Kekulize(molecule)
        self._init_from_rdkit_mol(
            molecule=molecule,
            functional_groups=functional_groups,
            identity_key=identity_key,
        )

    @classmethod
    def init_from_molecule(cls, molecule, functional_groups=None):
        """
        Initialize from a :class:`.Molecule`.

        Parameters
        ----------
        molecule : :class:`.Molecule`
            The molecule to initialize from.

        functional_groups : :class:`iterable` of :class:`str`, optional
            The names of the functional group types which are to be
            added to :attr:`func_groups`. If ``None``, then no
            functional groups are added.

        Returns
        -------
        :class:`.BuildingBlock`
             The building block. It will have the same atoms, bonds and
             atomic positions as `mol`.

        """

        return cls.init_from_rdkit_mol(
            molecule=molecule.to_rdkit_mol(),
            functional_groups=functional_groups,
        )

    @classmethod
    def init_from_file(cls, path, functional_groups=None):
        """
        Initialize from a file.

        Parameters
        ----------
        path : :class:`str`
            The path to a molecular structure file. Supported file
            types are:

                #. ``.mol``, ``.sdf`` - MDL V3000 MOL file
                #. ``.pdb`` - PDB file

        functional_groups : :class:`iterable` of :class:`str`, optional
            The names of the functional group types which are to be
            added to :attr:`func_groups`. If ``None``, then no
            functional groups are added.

        Returns
        -------
        :class:`.BuildingBlock`
            The building block.

        Raises
        ------
        :class:`ValueError`
            If the file type cannot be used for initialization.

        """

        if os.path.exists(path):
            _, ext = os.path.splitext(path)

            if ext not in cls._init_funcs:
                raise ValueError(
                    f'Unable to initialize from "{ext}" files.'
                )
            # This remake needs to be here because molecules loaded
            # with rdkit often have issues, because rdkit tries to do
            # bits of structural analysis like stereocenters. remake
            # gets rid of all this problematic metadata.
            molecule = remake(cls._init_funcs[ext](path))

        return cls.init_from_rdkit_mol(
            molecule=molecule,
            functional_groups=functional_groups,
        )

    @classmethod
    def init_from_random_file(
        cls,
        file_glob,
        functional_groups=None,
        random_seed=None,
    ):
        """
        Initialize from a random file in `file_glob`.

        Parameters
        ----------
        file_glob : :class:`str`
            A glob specifying files, one of which is used to initialize
            a :class:`.BuildingBlock` at random.

        functional_groups : :class:`iterable` of :class:`str`, optional
            The names of the functional group types which are to be
            added to :attr:`func_groups`. If ``None``, then no
            functional groups are added.

        random_seed : :class:`int`, optional
            The random seed to use.

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
        generator = np.random.RandomState(random_seed)
        generator.shuffle(files)

        for path in files:
            try:
                return cls.init_from_file(
                    path=path,
                    functional_groups=functional_groups,
                )

            except Exception:
                logger.warning(
                    f'Could not initialize {cls.__name__} from {path}.'
                )
        raise RuntimeError(
            f'No files in "{file_glob}" could be initialized from.'
        )

    @classmethod
    def init_from_rdkit_mol(cls, molecule, functional_groups=None):
        """
        Initialize from an :mod:`rdkit` molecule.

        Parameters
        ----------
        molecule : :class:`rdkit.Mol`
            The molecule.

        functional_groups : :class:`iterable` of :class:`str`, optional
            The names of the functional group types which are to be
            added to :attr:`func_groups`. If ``None``, then no
            functional groups are added.

        Returns
        -------
        :class:`BuildingBlock`
            The molecule.

        """

        key = cls._get_identity_key_from_rdkit_mol(
            molecule=molecule,
            functional_groups=functional_groups,
        )

        bb = cls.__new__(cls)
        bb._init_from_rdkit_mol(
            molecule=molecule,
            functional_groups=functional_groups,
            identity_key=key,
        )

        return bb

    def _init_from_rdkit_mol(
        self,
        molecule,
        functional_groups,
        identity_key,
    ):
        """
        Initialize from an :mod:`rdkit` molecule.

        Parameters
        ----------
        molecule : :class:`rdkit.Mol`
            The molecule.

        functional_groups : :class:`iterable` of :class:`str`
            The names of the functional group types which are to be
            added to :attr:`func_groups`. If ``None`, then no
            functional groups are added.

        identity_key : :class:`tuple`
            The identity key of the molecule.

        Returns
        -------
        None : :class:`NoneType`

        """

        if functional_groups is None:
            functional_groups = ()

        atoms = tuple(
            Atom(a.GetIdx(), a.GetAtomicNum(), a.GetFormalCharge())
            for a in molecule.GetAtoms()
        )
        bonds = tuple(
            Bond(
                atoms[b.GetBeginAtomIdx()],
                atoms[b.GetEndAtomIdx()],
                b.GetBondTypeAsDouble()
            )
            for b in molecule.GetBonds()
        )
        position_matrix = molecule.GetConformer().GetPositions()

        _Molecule.__init__(atoms, bonds, position_matrix)
        self._identity_key = identity_key

        fg_makers = (fg_types[name] for name in functional_groups)
        self._func_groups = tuple(
            func_group
            for fg_maker in fg_makers
            for func_group in fg_maker.get_functional_groups(self)
        )

    @classmethod
    def _init_from_dict(cls, molecule_dict):
        """
        Initialize from a :class:`dict` representation.

        Parameters
        ----------
        molecule_dict : :class:`dict`
            A :class:`dict` representation of a molecule generated
            by :meth:`to_dict`.

        Returns
        -------
        :class:`BuildingBlock`
            The molecule described by `mol_dict`.

        """

        d = dict(molecule_dict)
        identity_key = eval(d.pop('identity_key'))

        d.pop('class')
        functional_groups = d.pop('func_groups')

        obj = cls.__new__(cls)
        obj._identity_key = identity_key
        obj._position_matrix = np.array(d.pop('position_matrix')).T
        # If the cache is not being used, make sure to update all the
        # atoms and attributes to those in the dict.
        obj._atoms = eval(d.pop('atoms'), vars(atoms))
        obj._bonds = eval(d.pop('bonds'), vars(bonds))
        for bond in obj._bonds:
            bond.atom1 = obj.atoms[bond.atom1]
            bond.atom2 = obj.atoms[bond.atom2]

        fg_makers = (fg_types[name] for name in functional_groups)
        obj._func_groups = tuple(
            func_group
            for fg_maker in fg_makers
            for func_group in fg_maker.get_functional_groups(obj)
        )
        for attr, val in d.items():
            setattr(obj, attr, eval(val))

        return obj

    def clone(self):
        """
        Return a clone.

        Returns
        -------
        :class:`.BuildingBlock`
            The clone.

        """

        clone = super().clone()
        atom_map = {
            original: clone
            for original, clone in zip(self._atoms, clone._atoms)
        }
        clone._func_groups = tuple(
            fg.clone(atom_map) for fg in self._func_groups
        )
        return clone

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
            yield from self.func_groups[fg_id].get_bonder_ids()

    def get_bonder_centroids(self, fg_ids=None):
        """
        Yield the centroids of bonder atoms.

        A bonder centroid is the centroid of all bonder atoms in a
        particular functional group.

        Parameters
        ----------
        fg_ids : :class:`iterable` of :class:`int`
            The ids of functional groups to be used. The bonder
            centroids will be yielded in this order.
            If ``None`` then all functional groups are used and
            centroids are yielded in ascending order of functional
            group id.

        Yields
        ------
        :class:`numpy.ndarray`
            The centroid of a functional group.

        """

        if fg_ids is None:
            fg_ids = range(len(self.func_groups))

        for fg_id in fg_ids:
            yield self.get_centroid(
                atom_ids=self.func_groups[fg_id].get_bonder_ids()
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
            The ids of functional groups used to construct the plane.
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
            atom_ids=self.func_groups[fg_ids[0]].get_bonder_ids()
        )
        normal = self.get_bonder_plane_normal(fg_ids=fg_ids)
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
            The ids of functional groups used to construct the plane.
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
        :class:`ValueError`
            If there are not at least 3 functional groups, which is
            necessary to define a plane.

        """

        if fg_ids is None:
            fg_ids = range(len(self.func_groups))
        else:
            # The iterable is used mutliple times.
            fg_ids = list(fg_ids)

        if len(fg_ids) < 3:
            raise ValueError(
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
        if (
            # vector_angle is NaN if cc_vector is [0, 0, 0].
            not np.allclose(cc_vector, [0, 0, 0], atol=1e-5)
            and vector_angle(normal, cc_vector) > np.pi/2
        ):
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
                atom_ids=self.func_groups[i].get_bonder_ids()
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
                atom_ids=self.func_groups[i].get_bonder_ids()
            ))
            for i in fg_ids
        )
        pairs = it.combinations(iterable=centroids, r=2)
        for (id1, c1), (id2, c2) in pairs:
            yield id2, id1, c1-c2

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
            The vector running from the centroid of the bonder atoms to
            the molecular centroid.

        """

        if fg_ids is None:
            fg_ids = range(len(self.func_groups))
        else:
            # This iterable gets used more than once.
            fg_ids = list(fg_ids)

        bonder_centroid = self.get_centroid(
            atom_ids=self.get_bonder_ids(fg_ids=fg_ids)
        )
        return self.get_centroid() - bonder_centroid

    def to_dict(self, include_attrs=None, ignore_missing_attrs=False):
        """
        Return a :class:`dict` representation.

        Parameters
        ----------
        include_attrs : :class:`list` of :class:`str`, optional
            The names of additional attributes of the molecule to be
            added to the :class:`dict`. Each attribute is saved as a
            string using :func:`repr`.

        ignore_missing_attrs : :class:`bool`, optional
            If ``False`` and an attribute in `include_attrs` is not
            held by the :class:`BuildingBlock`, an error will be
            raised.

        Returns
        -------
        :class:`dict`
            A :class:`dict` which represents the molecule.

        """

        if include_attrs is None:
            include_attrs = []

        fgs = list(dedupe(fg.fg_type.name for fg in self.func_groups))

        bonds = []
        for bond in self.bonds:
            clone = bond.clone()
            clone.atom1 = clone.atom1.id
            clone.atom2 = clone.atom2.id
            bonds.append(clone)

        d = {
            'class': self.__class__.__name__,
            'func_groups': fgs,
            'position_matrix': self.get_position_matrix().tolist(),
            'atoms': repr(self.atoms),
            'bonds': repr(bonds),
            'identity_key': repr(self._identity_key)
        }

        if ignore_missing_attrs:
            d.update({
                attr: repr(getattr(self, attr))
                for attr in include_attrs
                if hasattr(self, attr)
            })
        else:
            d.update({
                attr: repr(getattr(self, attr))
                for attr in include_attrs
            })

        return d

    @staticmethod
    def _get_identity_key_from_rdkit_mol(mol, functional_groups):
        if functional_groups is None:
            functional_groups = ()
        functional_groups = sorted(functional_groups)

        # Don't modify the original molecule.
        mol = rdkit.Mol(mol)
        rdkit.SanitizeMol(mol)
        return (
            *functional_groups,
            rdkit.MolToSmiles(mol, canonical=True)
        )

    def __str__(self):
        smiles = rdkit.MolToSmiles(rdkit.RemoveHs(self.to_rdkit_mol()))
        func_groups = list(dedupe(
            fg.fg_type.name for fg in self.func_groups
        ))
        return f'{self.__class__.__name__}({smiles!r}, {func_groups})'

    def __repr__(self):
        return str(self)
