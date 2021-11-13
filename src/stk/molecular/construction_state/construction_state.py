"""
Construction State
==================

"""

from __future__ import annotations

import typing
from collections import abc
import numpy as np

from stk.utilities.typing import OneOrMany
from .graph_state import GraphState
from .molecule_state import MoleculeState
from ..edge import Edge
from ..edge_group import EdgeGroup
from ..vertex import Vertex
from ..placement_result import PlacementResult
from ..building_block import BuildingBlock
from ..atom import Atom
from ..atom_info import AtomInfo
from ..bond import Bond
from ..bond_info import BondInfo
from ..functional_groups import FunctionalGroup
from ..reactions import Reaction, ReactionResult


__all__ = (
    'ConstructionState',
)

_T = typing.TypeVar('_T', bound='ConstructionState')
#: A type variable matching any :class:`.Vertex`.
VertexT = typing.TypeVar('VertexT', bound=Vertex)

_LatticeConstants = tuple[np.ndarray, np.ndarray, np.ndarray]


class ConstructionState(typing.Generic[VertexT]):
    """
    The state of the molecule and topology graph under construction.

    """

    def __init__(
        self,
        building_block_vertices:
            dict[BuildingBlock, tuple[VertexT, ...]],
        edges: tuple[Edge, ...],
        lattice_constants: typing.Optional[_LatticeConstants] = None,
    ) -> None:
        """
        Initialize a :class:`.ConstructionState` instance.

        Parameters:

            building_block_vertices:
                Maps each :class:`.BuildingBlock` to be placed, to a
                :class:`tuple` of :class:`.Vertex` instances, on which
                it should be placed.

            edges:
                The edges of the topology graph.

            lattice_constants:
                A :class:`numpy.ndarray` for each lattice constant.
                Can be ``None`` if the topology graph is
                not periodic.

        """

        self._graph_state = GraphState(
            building_block_vertices=building_block_vertices,
            edges=edges,
            lattice_constants=lattice_constants,
        )
        self._molecule_state = MoleculeState()

    def _clone(self: _T) -> _T:
        clone = self.__class__.__new__(self.__class__)
        clone._graph_state = self._graph_state
        clone._molecule_state = self._molecule_state
        return clone

    def clone(self) -> ConstructionState[VertexT]:
        """
        Return a clone.

        Returns:

            The clone.

        """
        return self._clone()

    def _with_placement_results(
        self: _T,
        vertices: abc.Collection[VertexT],
        edges: abc.Iterable[abc.Sequence[Edge]],
        building_blocks: abc.Iterable[BuildingBlock],
        results: abc.Iterable[PlacementResult],
    ) -> _T:
        """
        Modify the instance.

        """

        self._molecule_state = (
            self._molecule_state.with_placement_results(
                vertices=vertices,
                edges=edges,
                building_blocks=building_blocks,
                results=results,
            )
        )
        return self

    def with_placement_results(
        self,
        vertices: abc.Collection[VertexT],
        edges: abc.Iterable[abc.Sequence[Edge]],
        building_blocks: abc.Iterable[BuildingBlock],
        results: abc.Iterable[PlacementResult],
    ) -> ConstructionState[VertexT]:
        """
        Return a clone holding the placement results.

        Parameters:

            vertices:
                The vertices used for placement.

            edges:
                For each vertex in `vertices`, the :class:`.Edge`
                instances connected to it.

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

    def get_lattice_constants(
        self,
    ) -> typing.Optional[_LatticeConstants]:
        """
        Get the lattice constants of the state.

        Returns:

            The lattice constants, or ``None`` if the topology graph
            is not periodic.

        """

        return self._graph_state.get_lattice_constants()

    def get_building_block(self, vertex_id: int) -> BuildingBlock:
        """
        Get the building block to be placed on a given vertex.

        Parameters:

            vertex_id:
                The id of the vertex, on which the building block is to
                be placed.

        Returns:

            The building block.

        """

        return self._graph_state.get_building_block(vertex_id)

    def get_vertices(
        self,
        vertex_ids: typing.Optional[OneOrMany[int]] = None,
    ) -> abc.Iterator[VertexT]:
        """
        Get some vertices.

        Parameters:

            vertex_ids:
                The ids of vertices to yield. If ``None``, all vertices
                will be yielded. Can be a single :class:`int` if a

        Yields:

            A vertex.

        """

        yield from self._graph_state.get_vertices(vertex_ids)

    def get_num_vertices(self) -> int:
        """
        Get the number of vertices in the topology graph.

        Returns:

            The number of vertices in the topology graph.

        """

        return self._graph_state.get_num_vertices()

    def get_edge(self, edge_id: int) -> Edge:
        """
        Get an edge.

        Parameters:

            edge_id:
                The id of an edge.

        Returns:

            An edge.

        """

        return self._graph_state.get_edge(edge_id)

    def get_num_edges(self) -> int:
        """
        Get the number of edges in the topology graph.

        Returns:

            The number of edges.

        """

        return self._graph_state.get_num_edges()

    def get_edges(self, vertex_id: int) -> abc.Sequence[Edge]:
        """
        Get the edges connected to a vertex.

        Parameters:

            vertex_id:
                The id of a vertex.

        Returns:

            The connected edges.

        """

        return self._graph_state.get_edges(vertex_id)

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

        yield from (
            self._molecule_state.get_edge_group_functional_groups(
                edge_group=edge_group,
            )
        )

    def _with_reaction_results(
        self: _T,
        reactions: abc.Collection[Reaction],
        results: abc.Iterable[ReactionResult],
    ) -> _T:
        """
        Modify the instance.

        """

        self._molecule_state = (
            self._molecule_state.with_reaction_results(
                reactions=reactions,
                results=results,
            )
        )
        return self

    def with_reaction_results(
        self,
        reactions: abc.Collection[Reaction],
        results: abc.Iterable[ReactionResult],
    ) -> ConstructionState[VertexT]:
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

    def _with_lattice_constants(
        self: _T,
        lattice_constants: typing.Optional[_LatticeConstants],
    ) -> _T:
        """
        Modify the instance.

        """

        self._graph_state = (
            self._graph_state.with_lattice_constants(
                lattice_constants=lattice_constants,
            )
        )
        return self

    def with_lattice_constants(
        self,
        lattice_constants: typing.Optional[_LatticeConstants],
    ) -> ConstructionState[VertexT]:
        """
        Return a clone holding the `lattice_constants`.

        Parameters:

            lattice_constants:
                The lattice constants of the clone. Requires 3 arrays,
                each of length 3, or ``None`` if the topology
                graph is not periodic.

        Returns:

            The clone holding the new lattice constants.

        """

        return self.clone()._with_lattice_constants(lattice_constants)

    def _with_position_matrix(
        self: _T,
        position_matrix: np.ndarray,
    ) -> _T:
        """
        Modify the instance.

        """

        self._molecule_state = (
            self._molecule_state.with_position_matrix(
                position_matrix=position_matrix,
            )
        )
        return self

    def with_position_matrix(
        self,
        position_matrix: np.ndarray,
    ) -> ConstructionState[VertexT]:
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

    def _with_vertices(
        self: _T,
        vertices: abc.Iterable[VertexT],
    ) -> _T:
        """
        Modify the instance.

        """

        self._graph_state = self._graph_state.with_vertices(vertices)
        return self

    def with_vertices(
        self,
        vertices: abc.Iterable[VertexT],
    ) -> ConstructionState[VertexT]:
        """
        Returns a clone holding `vertices`.

        Parameters:

            vertices:
                The vertices the clone should hold.

        Returns:

            The clone.

        """

        return self.clone()._with_vertices(vertices)

    def get_position_matrix(self) -> np.ndarray:
        """
        Get the position matrix of the molecule being constructed.

        Returns:

            The position matrix.

        """

        return self._molecule_state.get_position_matrix()

    def get_atoms(self) -> abc.Iterator[Atom]:
        """
        Yield the atoms of the molecule being constructed.

        Yields:

            An atom of the molecule being constructed.

        """

        yield from self._molecule_state.get_atoms()

    def get_bonds(self) -> abc.Iterator[Bond]:
        """
        Yield the bonds of the molecule being constructed.

        Yields:

            A bond of the molecule being constructed.

        """

        yield from self._molecule_state.get_bonds()

    def get_atom_infos(self) -> abc.Iterator[AtomInfo]:
        """
        Yield the atom infos of the molecule being constructed.

        Yields:

            An atom info of the molecule being constructed.

        """

        yield from self._molecule_state.get_atom_infos()

    def get_bond_infos(self) -> abc.Iterator[BondInfo]:
        """
        Yield the bond infos of the molecule being constructed.

        Yields:

            The bond info of the molecule being constructed.

        """

        yield from self._molecule_state.get_bond_infos()

    def get_num_building_block(
        self,
        building_block: BuildingBlock,
    ) -> int:
        """
        Get the number of times `building_block` is present.

        Parameters:

            building_block:
                The building block whose frequency in the topology
                graph is desired.

        Returns:

            The number of times `building_block` is present in the
            topology graph.

        """

        return self._graph_state.get_num_building_block(building_block)

    def get_building_blocks(self) -> abc.Iterator[BuildingBlock]:
        """
        Yield the building blocks.

        Building blocks are yielded in an order based on their
        position in the topology graph. For two equivalent
        topology graphs, but with different building blocks,
        equivalently positioned building blocks will be yielded at the
        same time.

        Yields:

            A building block of the topology graph.

        """

        yield from self._graph_state.get_building_blocks()
