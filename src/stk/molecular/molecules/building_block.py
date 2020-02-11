"""
Building Block
==============

"""


import logging
import os
import numpy as np
import rdkit.Chem.AllChem as rdkit
from functools import partial

from .. import atoms, functional_groups
from ..functional_groups import FunctionalGroup
from ..atoms import Atom
from ..bond import Bond
from .molecule import Molecule_
from ...utilities import remake


logger = logging.getLogger(__name__)


class BuildingBlock(Molecule_):
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
        rdkit.Kekulize(molecule)
        self._init_from_rdkit_mol(molecule, functional_groups)

    @classmethod
    def init_from_molecule(cls, molecule, functional_groups=None):
        """
        Initialize from a :class:`.Molecule`.

        Parameters
        ----------
        molecule : :class:`.Molecule`
            The molecule to initialize from.

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
    def init_from_rdkit_mol(cls, molecule, functional_groups=None):
        """
        Initialize from an :mod:`rdkit` molecule.

        Parameters
        ----------
        molecule : :class:`rdkit.Mol`
            The molecule.

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

        Returns
        -------
        :class:`BuildingBlock`
            The molecule.

        """

        bb = cls.__new__(cls)
        bb._init_from_rdkit_mol(
            molecule=molecule,
            functional_groups=functional_groups,
        )

        return bb

    def _init_from_rdkit_mol(
        self,
        molecule,
        functional_groups,
    ):
        """
        Initialize from an :mod:`rdkit` molecule.

        Parameters
        ----------
        molecule : :class:`rdkit.Mol`
            The molecule.

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
                atom1=atoms[b.GetBeginAtomIdx()],
                atom2=atoms[b.GetEndAtomIdx()],
                order=b.GetBondTypeAsDouble(),
            )
            for b in molecule.GetBonds()
        )
        position_matrix = molecule.GetConformer().GetPositions()

        super().__init__(atoms, bonds, position_matrix)
        self._with_functional_groups(self._extract_functional_groups(
            functional_groups=functional_groups,
        ))

    def _extract_functional_groups(self, functional_groups):
        """
        Yield functional groups.

        The input can be a mixture of :class:`.FunctionalGroup` and
        :class:`.FunctionalGroupFactory`. The output yields
        :class:`.FunctionalGroup` instances only. Either those
        held directly in `functional_groups` or created by the
        factories in `functional_groups`.

        Parameters
        ----------
        functional_groups : :class:`iterable`
            Can be an :class:`iterable` of both
            :class:`.FunctionalGroup` and
            :class:`.FunctionalGroupFactory`.

        Yields
        ------
        :class:`.FunctionalGroup`
            A functional group from `functional_groups`, or created
            by a factory in `functional_groups`.

        """

        for functional_group in functional_groups:
            if isinstance(functional_group, FunctionalGroup):
                yield functional_group
            else:
                # Else it's a factory.
                yield from functional_group.get_functional_groups(self)

    def _with_functional_groups(self, functional_groups):
        """
        Modify the molecule.

        """

        atom_map = {a.get_id(): a for a in self._atoms}
        self._functional_groups = tuple(
            fg.with_atoms(atom_map) for fg in functional_groups
        )
        return self

    def with_functional_groups(self, functional_groups):
        """
        Return a clone with specific functional groups.

        Parameters
        ----------
        functional_groups : :class:`iterable`
            :class:`.FunctionalGroup` instances which the clone
            should have.

        Returns
        -------
        :class:`.BuildingBlock`
            The clone.

        """

        return self.clone()._with_functional_groups(functional_groups)

    def get_num_functional_groups(self):
        """
        Return the number of functional groups.

        Returns
        -------
        :class:`int`
            The number of functional groups in the building block.

        """

        return len(self._functional_groups)

    def get_functional_groups(self, fg_ids=None):
        """
        Yield the functional groups, ordered by id.

        Parameters
        ----------
        fg_ids : :class:`iterable` of :class:`int`, optional
            The ids of functional groups yielded. If ``None``, then
            all functional groups are yielded. Can be a single
            :class:`int`, if a single functional group is
            desired.

        Yields
        ------
        :class:`.FunctionalGroup`
            A functional group of the building block.

        """

        if fg_ids is None:
            fg_ids = range(len(self._functional_groups))
        elif isinstance(fg_ids, int):
            fg_ids = (fg_ids, )

        for fg_id in fg_ids:
            yield self._functional_groups[fg_id]

    @classmethod
    def init_from_dict(cls, molecule_dict):
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
            The molecule described by `molecule_dict`.

        """

        obj = super().init_from_dict(molecule_dict)
        obj = cls.__new__(cls)
        obj._position_matrix = np.array(
            molecule_dict['position_matrix']
        ).T
        obj._atoms = eval(molecule_dict['atoms'], vars(atoms))
        obj._bonds = [
            Bond(
                atom1=obj._atoms[bond_dict['atom1_id']],
                atom2=obj._atoms[bond_dict['atom2_id']],
                order=bond_dict['order'],
                periodicity=tuple(bond_dict['periodicity']),
            )
            for bond_dict in molecule_dict['bonds']
        ]

        obj._functional_groups = []
        globals_ = vars(functional_groups)
        globals_.update(vars(atoms))
        obj._with_functional_groups(
            functional_groups=eval(
                molecule_dict['functional_groups'],
                globals_
            )
        )
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
        clone._functional_groups = self._functional_groups
        return clone

    def get_bonder_ids(self, fg_ids=None):
        """
        Yield ids of bonder atoms.

        Parameters
        ----------
        fg_ids : :class:`iterable` of :class:`int`
            The ids of functional groups whose bonder atoms should be
            yielded. Can be a single :class:`int`, if a single
            functional group should be used, or  ``None``, if all
            functional groups should be used.

        Yields
        ------
        :class:`int`
            The id of a bonder atom.

        """

        if fg_ids is None:
            fg_ids = range(len(self._functional_groups))
        elif isinstance(fg_ids, int):
            fg_ids = (fg_ids, )

        for fg_id in fg_ids:
            yield from self._functional_groups[fg_id].get_bonder_ids()

    def to_dict(self):
        """
        Return a :class:`dict` representation.

        Returns
        -------
        :class:`dict`
            A :class:`dict` which represents the molecule.

        """

        return {
            'class': self.__class__.__name__,
            'functional_groups': repr(self._functional_groups),
            'position_matrix': self.get_position_matrix().tolist(),
            'atoms': repr(self._atoms),
            'bonds': [b.to_dict() for b in self._bonds],
        }

    def __str__(self):
        if self._functional_groups:
            fg_repr = f', {self._functional_groups!r}'
        else:
            fg_repr = ''

        smiles = rdkit.MolToSmiles(rdkit.RemoveHs(self.to_rdkit_mol()))
        return f'{self.__class__.__name__}({smiles!r}{fg_repr})'

    def __repr__(self):
        return str(self)
