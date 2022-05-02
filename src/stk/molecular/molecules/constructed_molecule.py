"""
Constructed Molecule
====================

"""

from __future__ import annotations

import logging
import typing

import numpy as np
import rdkit.Chem.AllChem as rdkit

from ...utilities import OneOrMany
from ..atoms import Atom, AtomInfo
from ..bonds import Bond, BondInfo
from ..topology_graphs import TopologyGraph
from ..topology_graphs.topology_graph.construction_result import (
    ConstructionResult,
)
from .molecule import Molecule
from .utilities import get_bond_info_atom_ids, sort_bond_atoms_by_id

logger = logging.getLogger(__name__)


class ConstructedMolecule(Molecule):
    """
    Represents constructed molecules.

    Examples:

        *Initialization*

        A :class:`.ConstructedMolecule` is initialized from a
        :class:`.TopologyGraph`, which is typically initialized from
        some :class:`.BuildingBlock` instances.

        .. testcode:: initialization

            import stk

            bb1 = stk.BuildingBlock(
                smiles='NCCCN',
                functional_groups=[stk.PrimaryAminoFactory()],
            )
            bb2 = stk.BuildingBlock(
                smiles='O=CC(C=O)CC=O',
                functional_groups=[stk.AldehydeFactory()],
            )
            tetrahedron = stk.cage.FourPlusSix((bb1, bb2))
            cage = stk.ConstructedMolecule(tetrahedron)

        *Hierarchical Construction*

        A :class:`ConstructedMolecule` may be used to construct other
        :class:`ConstructedMolecule` instances, though you will
        probably have to convert it to a :class:`.BuildingBlock` first

        .. testcode:: hierarchical-construction

            import stk

            bb1 = stk.BuildingBlock(
                smiles='NCCCN',
                functional_groups=[stk.PrimaryAminoFactory()],
            )
            bb2 = stk.BuildingBlock(
                smiles='O=CC(C=O)CC=O',
                functional_groups=[stk.AldehydeFactory()],
            )
            tetrahedron = stk.cage.FourPlusSix((bb1, bb2))
            cage = stk.ConstructedMolecule(tetrahedron)

            benzene = stk.host_guest.Guest(
                building_block=stk.BuildingBlock('c1ccccc1'),
            )
            cage_complex = stk.host_guest.Complex(
                host=stk.BuildingBlock.init_from_molecule(cage),
                guests=benzene,
            )
            cage_complex = stk.ConstructedMolecule(cage_complex)

        Obviously, the initialization of the
        :class:`.ConstructedMolecule` depends mostly on the specifics
        of the :class:`.TopologyGraph` used, and the documentation of
        those classes should be examined for more examples.

    """

    _atom_infos:  tuple[AtomInfo, ...]
    _bond_infos: tuple[BondInfo, ...]
    _num_building_blocks: dict[Molecule, int]

    def __init__(
        self,
        topology_graph: TopologyGraph,
    ) -> None:
        """
        Initialize a :class:`.ConstructedMolecule`.

        Parameters:

            topology_graph:
                The topology graph of the constructed molecule.

        """

        self._init_from_construction_result(
            obj=self,
            construction_result=topology_graph.construct(),
        )

    @classmethod
    def init(
        cls,
        atoms: tuple[Atom, ...],
        bonds: tuple[Bond, ...],
        position_matrix: np.ndarray,
        atom_infos: tuple[AtomInfo, ...],
        bond_infos: tuple[BondInfo, ...],
        num_building_blocks: dict[Molecule, int],
    ) -> ConstructedMolecule:
        """
        Initialize a :class:`.ConstructedMolecule` from its components.

        Parameters:

            atoms:
                The atoms of the molecule.

            bond:
                The bonds of the molecule.

            position_matrix:
                A ``(n, 3)`` position matrix of the molecule.

            atom_infos:
                The atom infos of the molecule.

            bond_infos:
                The bond infos of the molecule.

            num_building_blocks:
                Maps each building block of the constructed molecule to
                the number of times it is present in it.

        Returns:

            The constructed molecule.

        """

        molecule = cls.__new__(cls)
        Molecule.__init__(molecule, atoms, bonds, position_matrix)
        molecule._atom_infos = atom_infos
        molecule._bond_infos = bond_infos
        molecule._num_building_blocks = dict(num_building_blocks)
        return molecule

    @classmethod
    def init_from_construction_result(
        cls,
        construction_result: ConstructionResult,
    ) -> ConstructedMolecule:
        """
        Initialize a :class:`.ConstructedMolecule`.

        Parameters:

            construction_result:
                The result of a construction, from which the
                :class:`.ConstructedMolecule` should be initialized.

        Returns:

            The constructed molecule.

        """

        return cls._init_from_construction_result(
            obj=cls.__new__(cls),
            construction_result=construction_result,
        )

    @staticmethod
    def _init_from_construction_result(
        obj: ConstructedMolecule,
        construction_result: ConstructionResult,
    ) -> ConstructedMolecule:
        """
        Initialize a :class:`.ConstructedMolecule`.

        This modifies `obj`.

        Parameters:

            obj:
                The constructed molecule to initialize.

            construction_result:
                The result of a construction, from which the
                :class:`.ConstructedMolecule` should be initialized.

        Returns:

            The `obj` instance.

        """

        super(ConstructedMolecule, obj).__init__(
            atoms=construction_result.get_atoms(),
            bonds=construction_result.get_bonds(),
            position_matrix=construction_result.get_position_matrix(),
        )
        obj._atom_infos = construction_result.get_atom_infos()
        obj._bond_infos = construction_result.get_bond_infos()
        obj._num_building_blocks = {
            building_block: construction_result.get_num_building_block(
                building_block=building_block,
            )
            for building_block
            in construction_result.get_building_blocks()
        }
        return obj

    def clone(self) -> ConstructedMolecule:
        clone = self._clone()
        clone._atom_infos = self._atom_infos
        clone._bond_infos = self._bond_infos
        clone._num_building_blocks = dict(self._num_building_blocks)
        return clone

    def get_building_blocks(self) -> typing.Iterator[Molecule]:
        """
        Yield the building blocks of the constructed molecule.

        Building blocks are yielded in an order based on their
        position in the constructed molecule. For two topologically
        equivalent constructed molecules, but with different building
        blocks, equivalently positioned building blocks will be
        yielded at the same time.

        Yields:

            A building block of the constructed molecule.

        """

        yield from self._num_building_blocks

    def get_num_building_block(
        self,
        building_block: Molecule,
    ) -> int:
        """
        Get the number of times `building_block` is present.

        Parameters:

            building_block:
                The building block whose frequency in the constructed
                molecule is desired.

        Returns:

            The number of times `building_block` was used in the
            construction of the constructed molecule.

        """

        return self._num_building_blocks[building_block]

    def get_atom_infos(
        self,
        atom_ids: typing.Optional[OneOrMany[int]] = None,
    ) -> typing.Iterator[AtomInfo]:
        """
        Yield data about atoms in the molecule.

        Parameters:

            atom_ids:
                The ids of atoms whose data is desired. If ``None``,
                data on all atoms will be yielded. Can be a single
                :class:`int`, if data on a single atom is desired.

        Yields:

            Data about an atom.

        """

        if atom_ids is None:
            atom_ids = range(len(self._atoms))
        elif isinstance(atom_ids, int):
            atom_ids = (atom_ids, )

        for atom_id in atom_ids:
            yield self._atom_infos[atom_id]

    def get_bond_infos(self) -> typing.Iterator[BondInfo]:
        """
        Yield data about bonds in the molecule.

        Yields:

            Data about a bond.

        """

        yield from self._bond_infos

    def with_canonical_atom_ordering(self) -> ConstructedMolecule:
        return self.clone()._with_canonical_atom_ordering()

    def _with_canonical_atom_ordering(self) -> ConstructedMolecule:
        # Make all building blocks canonically ordered too.
        building_blocks = {
            building_block:
                building_block.with_canonical_atom_ordering()

            for building_block in self._num_building_blocks
        }

        # Cache these mappings for later, to avoid unnecessary
        # re-computations of canonical ordering.
        canonical_map = {
            building_block: building_block.get_canonical_atom_ids()
            for building_block in self._num_building_blocks
        }

        self._num_building_blocks = {
            building_block: num
            for building_block, num
            in zip(
                building_blocks.values(),
                self._num_building_blocks.values(),
            )
        }

        ordering = rdkit.CanonicalRankAtoms(self.to_rdkit_mol())
        id_map = {
            new_id: atom.get_id()
            for new_id, atom in zip(ordering, self._atoms)
        }
        super()._with_canonical_atom_ordering()
        atom_map = {
            old_id: self._atoms[new_id]
            for old_id, new_id in enumerate(ordering)
        }
        old_atom_infos = self._atom_infos

        def get_atom_info(atom: Atom) -> AtomInfo:

            old_atom_info = old_atom_infos[id_map[atom.get_id()]]
            old_building_block = old_atom_info.get_building_block()

            if old_building_block is None:
                return AtomInfo(
                    atom=atom,
                    building_block_atom=None,
                    building_block=None,
                    building_block_id=None,
                )

            old_building_block_atom = (
                old_atom_info.get_building_block_atom()
            )

            canonical_building_block_atom_id = canonical_map[
                old_building_block
            ][old_building_block_atom.get_id()]

            canonical_building_block = building_blocks[
                old_building_block
            ]

            canonical_building_block_atom, = (
                canonical_building_block.get_atoms(
                    atom_ids=canonical_building_block_atom_id,
                )
            )

            return AtomInfo(
                atom=atom,
                building_block_atom=canonical_building_block_atom,
                building_block=canonical_building_block,
                building_block_id=(
                    old_atom_info.get_building_block_id()
                ),
            )

        def get_bond_info(info: BondInfo) -> BondInfo:
            building_block = info.get_building_block()
            return BondInfo(
                bond=sort_bond_atoms_by_id(
                    info.get_bond().with_atoms(atom_map)
                ),
                building_block=(
                    building_block
                    if building_block is None
                    else building_blocks[building_block]
                ),
                building_block_id=info.get_building_block_id(),
            )

        self._atom_infos = tuple(map(get_atom_info, self._atoms))
        self._bond_infos = tuple(sorted(
            map(get_bond_info, self._bond_infos),
            key=get_bond_info_atom_ids,
        ))
        return self

    def with_centroid(
        self,
        position: np.ndarray,
        atom_ids: typing.Optional[OneOrMany[int]] = None,
    ) -> ConstructedMolecule:

        return self.clone()._with_centroid(position, atom_ids)

    def with_displacement(
        self,
        displacement: np.ndarray,
    ) -> ConstructedMolecule:

        return self.clone()._with_displacement(displacement)

    def with_position_matrix(
        self,
        position_matrix: np.ndarray,
    ) -> ConstructedMolecule:

        return self.clone()._with_position_matrix(position_matrix)

    def with_rotation_about_axis(
        self,
        angle: float,
        axis: np.ndarray,
        origin: np.ndarray,
    ) -> ConstructedMolecule:

        return self.clone()._with_rotation_about_axis(
            angle=angle,
            axis=axis,
            origin=origin,
        )

    def with_rotation_between_vectors(
        self,
        start: np.ndarray,
        target: np.ndarray,
        origin: np.ndarray,
    ) -> ConstructedMolecule:

        return self.clone()._with_rotation_between_vectors(
            start=start,
            target=target,
            origin=origin,
        )

    def with_rotation_to_minimize_angle(
        self,
        start: np.ndarray,
        target: np.ndarray,
        axis: np.ndarray,
        origin: np.ndarray,
    ) -> ConstructedMolecule:

        return self.clone()._with_rotation_to_minimize_angle(
            start=start,
            target=target,
            axis=axis,
            origin=origin,
        )

    def with_structure_from_file(
        self,
        path: str,
        extension: typing.Optional[str] = None,
    ) -> ConstructedMolecule:

        return typing.cast(
            ConstructedMolecule,
            super().with_structure_from_file(path, extension),
        )

    def write(
        self,
        path: str,
        atom_ids: typing.Optional[OneOrMany[int]] = None,
    ) -> ConstructedMolecule:

        return typing.cast(
            ConstructedMolecule,
            super().write(path, atom_ids),
        )
