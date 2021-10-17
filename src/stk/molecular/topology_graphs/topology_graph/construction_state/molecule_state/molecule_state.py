"""
Molecule State
==============

"""

from __future__ import annotations

__all__ = (
    'MoleculeState',
)

import numpy as np
import typing
from collections import abc

from .reactions_summary import ReactionsSummary
from .deletions_summary import DeletionsSummary
from .placements_summary import PlacementsSummary
from ...vertex import Vertex
from ...edge import Edge
from ...edge_group import EdgeGroup
from ...placement_result import PlacementResult
from .....building_block import BuildingBlock
from .....functional_groups import FunctionalGroup
from .....atom import Atom
from .....atom_info import AtomInfo
from .....bond import Bond
from .....bond_info import BondInfo
from .....reactions import Reaction, ReactionResult


_T = typing.TypeVar('_T', bound='MoleculeState')


class MoleculeState:
    """
    Represents the state of a molecule under construction.

    """

    _atoms: list[Atom]
    _atom_infos: list[AtomInfo]
    _bonds: list[Bond]
    _bond_infos: list[BondInfo]
    _edge_functional_groups: dict[int, list[FunctionalGroup]]
    _num_placements: int

    __slots__ = [
        '_position_matrix',
        '_atoms',
        '_atom_infos',
        '_bonds',
        '_bond_infos',
        '_edge_functional_groups',
        '_num_placements',
    ]

    def __init__(self) -> None:
        """
        Initialize a :class:`.MoleculeState` instance.

        """

        self._position_matrix = np.empty((0, 3), dtype=np.float64)
        self._atoms = []
        self._atom_infos = []
        self._bonds = []
        self._bond_infos = []
        self._edge_functional_groups = {}
        self._num_placements = 0

    def clone(self) -> MoleculeState:
        """
        Return a clone.

        Returns:

            The clone.

        """

        clone = self.__class__.__new__(self.__class__)
        clone._position_matrix = np.array(self._position_matrix)
        clone._atoms = list(self._atoms)
        clone._atom_infos = list(self._atom_infos)
        clone._bonds = list(self._bonds)
        clone._bond_infos = list(self._bond_infos)
        clone._num_placements = self._num_placements
        clone._edge_functional_groups = {
            edge_id: list(functional_groups)
            for edge_id, functional_groups
            in self._edge_functional_groups.items()
        }
        return clone

    def _with_placement_results(
        self: _T,
        vertices: tuple[Vertex, ...],
        edges: abc.Iterable[abc.Sequence[Edge]],
        building_blocks: tuple[BuildingBlock, ...],
        results: abc.Iterable[PlacementResult],
    ) -> _T:
        """
        Modify the instance.

        """

        summary = PlacementsSummary(
            building_blocks=building_blocks,
            placement_results=results,
            num_atoms=len(self._atoms),
            num_previous_placements=self._num_placements,
        )
        self._num_placements += len(vertices)
        self._atoms.extend(summary.get_atoms())
        self._atom_infos.extend(summary.get_atom_infos())
        self._bonds.extend(summary.get_bonds())
        self._bond_infos.extend(summary.get_bond_infos())
        self._position_matrix = np.concatenate([
            self._position_matrix,
            summary.get_position_matrix(),
        ])
        for edge_id, functional_groups in (
            summary.get_edge_functional_groups()
        ):
            self._edge_functional_groups[edge_id] = (
                self._edge_functional_groups.get(edge_id, [])
            )
            self._edge_functional_groups[edge_id].extend(
                functional_groups
            )
        return self

    def with_placement_results(
        self,
        vertices: tuple[Vertex, ...],
        edges: abc.Iterable[abc.Sequence[Edge]],
        building_blocks: tuple[BuildingBlock, ...],
        results: abc.Iterable[PlacementResult],
    ) -> MoleculeState:
        """
        Return a clone holding the placement results.

        Parameters:

            vertices:
                The vertices used for placement.

            edges:
                For each vertex in `vertices`, a :class:`tuple` of
                :class:`.Edge` instances connected to it.

            building_blocks:
                For each vertex in `vertices`, the building block
                placed on it.

            results:
                For every vertex in `vertices`, the result of the
                placement.

        Returns:

            The clone holding the placement results.

        """

        return self.clone()._with_placement_results(
            vertices=vertices,
            edges=edges,
            building_blocks=building_blocks,
            results=results,
        )

    def get_position_matrix(self) -> np.ndarray:
        """
        Get the position matrix of the molecule.

        Returns:

            The position matrix.

        """

        return np.array(self._position_matrix)

    def get_atoms(self) -> abc.Iterator[Atom]:
        """
        Yield the atoms of the molecule.

        Yields:

            An atom of the molecule.

        """

        yield from self._atoms

    def get_bonds(self) -> abc.Iterator[Bond]:
        """
        Yield the bonds of the molecule.

        Yields
        ------
        :class:`.Bond`
            A bond of the molecule.

        """

        yield from self._bonds

    def get_atom_infos(self) -> abc.Iterator[AtomInfo]:
        """
        Yield the atom infos of the molecule.

        Yields:

            An atom info of the molecule.

        """

        yield from self._atom_infos

    def get_bond_infos(self) -> abc.Iterator[BondInfo]:
        """
        Yield the bond infos of the molecule.

        Yields:

            The bond info of the molecule.

        """

        yield from self._bond_infos

    def get_edge_group_functional_groups(
        self,
        edge_group: EdgeGroup,
    ) -> abc.Iterator[FunctionalGroup]:
        """
        Yield the functional groups associated with `edge_group`.

        Parameters:

            edge_group:
                The edge group, whose functional groups are desired.

        Yields:

            A functional group which belongs to `edge_group`.

        """

        for edge_id in edge_group.get_edge_ids():
            yield from self._edge_functional_groups[edge_id]

    def with_position_matrix(
        self,
        position_matrix: np.ndarray,
    ) -> MoleculeState:
        """
        Return a clone holding the `position_matrix`.

        Parameters:

            position_matrix:
                The position matrix of the clone. The shape of the
                matrix is ``(n, 3)``.

        Returns:

            The clone holding the new position matrix.

        """

        return self.clone()._with_position_matrix(position_matrix)

    def _with_position_matrix(
        self: _T,
        position_matrix: np.ndarray,
    ) -> _T:
        """
        Modify the instance.

        """

        self._position_matrix = np.array(position_matrix)
        return self

    def with_reaction_results(
        self,
        reactions: tuple[Reaction, ...],
        results: abc.Iterable[ReactionResult],
    ) -> MoleculeState:
        """
        Return a clone holding the reaction results.

        Parameters:

            reactions:
                The reactions.

            results:
                For each reaction in `reactions`, its result.

        Returns:

            The clone holding the reaction results.

        """

        return self.clone()._with_reaction_results(reactions, results)

    def _with_reaction_results(
        self: _T,
        reactions: tuple[Reaction, ...],
        results: abc.Iterable[ReactionResult],
    ) -> _T:
        """
        Modify the instance.

        """

        reactions_summary = ReactionsSummary(
            num_atoms=len(self._atoms),
            reaction_results=results,
        )
        self._with_reactions_summary(reactions_summary)
        self._with_deletions_summary(DeletionsSummary(
            atoms=self._atoms,
            atom_infos=self._atom_infos,
            bonds=self._bonds,
            bond_infos=self._bond_infos,
            position_matrix=self._position_matrix,
            deleted_atom_ids=reactions_summary.get_deleted_atom_ids(),
            deleted_bond_ids=reactions_summary.get_deleted_bond_ids(),
        ))
        return self

    def _with_reactions_summary(
        self,
        summary: ReactionsSummary,
    ) -> None:
        """
        Add the results held in `summary`.

        Parameters:

            summary:
                A summary of the reaction results.

        """

        self._atoms.extend(summary.get_atoms())
        self._atom_infos.extend(summary.get_atom_infos())
        self._bonds.extend(summary.get_bonds())
        self._bond_infos.extend(summary.get_bond_infos())

        positions = tuple(summary.get_positions())
        if positions:
            self._position_matrix = np.vstack([
                self._position_matrix,
                positions,
            ])

    def _with_deletions_summary(
        self: _T,
        summary: DeletionsSummary,
    ) -> _T:
        """
        Add the results held in `summary`.

        Parameters:

            summary:
                A summary of the atom deletion results.

        """

        self._atoms = list(summary.get_atoms())
        self._atom_infos = list(summary.get_atom_infos())
        self._bonds = list(summary.get_bonds())
        self._bond_infos = list(summary.get_bond_infos())
        self._position_matrix = np.array(list(summary.get_positions()))
        return self
