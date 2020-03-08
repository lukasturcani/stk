"""
Building Block
==============

"""


import logging
import os
import rdkit.Chem.AllChem as rdkit
from functools import partial

from ..functional_groups import FunctionalGroup
from ..atoms import Atom
from ..bond import Bond
from .molecule import Molecule
from ...utilities import remake


logger = logging.getLogger(__name__)


class BuildingBlock(Molecule):
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
            removeHs=False,
            proximityBonding=False,
        ),
    }

    def __init__(self, smiles, functional_groups=(), random_seed=4):
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
    def init_from_molecule(cls, molecule, functional_groups=()):
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

        return cls.init(
            atoms=tuple(molecule.get_atoms()),
            bonds=tuple(molecule.get_bonds()),
            position_matrix=molecule.get_position_matrix(),
            functional_groups=functional_groups,
        )

    @classmethod
    def init(
        cls,
        atoms,
        bonds,
        position_matrix,
        functional_groups,
    ):
        building_block = cls.__new__(cls)
        Molecule.__init__(
            self=building_block,
            atoms=atoms,
            bonds=bonds,
            position_matrix=position_matrix,
        )
        building_block._with_functional_groups(functional_groups)
        return building_block

    @classmethod
    def init_from_file(cls, path, functional_groups=()):
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
    def init_from_rdkit_mol(cls, molecule, functional_groups=()):
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

    def get_placer_ids(self):
        """
        Yield ids of atoms used for placing the building block.

        Yields
        ------
        :class:`int`
            The id of a placer atom.

        """

        for functional_group in self._functional_groups:
            yield from functional_group.get_placer_ids()

    def get_core_atom_ids(self):
        """
        Yield ids of atoms which form the core of the building block.

        This includes all atoms in the building block not part of a
        functional group, as well as any atoms in a functional group,
        specifically labelled as core atoms.

        Yields
        ------
        :class:`int`
            The id of a core atom.

        """

        functional_group_atom_ids = {
            atom_id
            for functional_group in self._functional_groups
            for atom_id in functional_group.get_atom_ids()
        }
        for atom in self._atoms:
            atom_id = atom.get_id()
            if atom_id not in functional_group_atom_ids:
                yield atom_id

        for functional_group in self._functional_groups:
            for atom_id in functional_group.get_core_atom_ids():
                yield atom_id

    def __str__(self):
        if self._functional_groups:
            fg_repr = f', {self._functional_groups!r}'
        else:
            fg_repr = ''

        smiles = rdkit.MolToSmiles(
            mol=rdkit.RemoveHs(self.to_rdkit_mol()),
            # Canonical smiles generation can be slow.
            canonical=False,
        )
        return f'{self.__class__.__name__}({smiles!r}{fg_repr})'

    def __repr__(self):
        return str(self)
