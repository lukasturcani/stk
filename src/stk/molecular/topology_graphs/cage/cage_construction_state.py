"""
Cage Construction State
=======================

"""

from __future__ import annotations

import typing
from collections import abc
import numpy as np

from .vertices import CageVertex
from ...building_block import BuildingBlock
from ...vertex import Vertex
from ...edge import Edge
from ...placement_result import PlacementResult
from ...construction_state import ConstructionState

__all__ = (
    'CageConstructionState',
)

_T = typing.TypeVar('_T', bound='CageConstructionState')


class CageConstructionState(ConstructionState[CageVertex]):
    """
    The construction state of a :class:`.Cage`.

    """

    _num_placement_stages_done: int
    _neighbor_positions: dict[int, list[np.ndarray]]

    def __init__(
        self,
        building_block_vertices:
            dict[BuildingBlock, tuple[CageVertex, ...]],
        edges: tuple[Edge, ...],
        num_placement_stages: int,
        vertex_degrees: dict[int, int],
    ) -> None:
        """
        Initialize a :class:`.CageConstructionState` instance.

        Parameters:

            building_block_vertices:
                Maps each :class:`.BuildingBlock` to be placed, to a
                :class:`tuple` of :class:`.Vertex` instances, on which
                it should be placed.

            edges:
                The edges of the topology graph.

            num_placement_stages:
                The number of placement stages.

            vertex_degrees:
                Maps the id of a vertex to the number of edges it is
                connected to.

        """

        super().__init__(
            building_block_vertices=building_block_vertices,
            edges=edges,
            lattice_constants=None,
        )
        self._num_placement_stages = num_placement_stages
        self._num_placement_stages_done = 0
        self._vertex_degrees = dict(vertex_degrees)
        self._neighbor_positions = {}

    def _with_placement_results(
        self: _T,
        vertices: tuple[CageVertex, ...],
        edges: abc.Iterable[abc.Sequence[Edge]],
        building_blocks: tuple[BuildingBlock, ...],
        results: abc.Iterable[PlacementResult],
    ) -> _T:
        # Need to iterate multiple times through results.
        results = tuple(results)
        super()._with_placement_results(
            vertices=vertices,
            edges=edges,
            building_blocks=building_blocks,
            results=results,
        )

        # No need to update vertex positions after the last stage.
        self._num_placement_stages_done += 1
        if (
            self._num_placement_stages_done
            >= self._num_placement_stages
        ):
            return self

        self._update_neighbor_positions(
            vertices=vertices,
            edges=edges,
            building_blocks=building_blocks,
            results=results,
        )
        return self

    def _update_neighbor_positions(
        self,
        vertices: abc.Iterable[CageVertex],
        edges: abc.Iterable[abc.Collection[Edge]],
        building_blocks: tuple[BuildingBlock, ...],
        results: abc.Iterable[PlacementResult],
    ) -> None:
        # Use literal docstring here to prevent linter errors stemming
        # from the use of "\" in the docstring.
        r"""
        Put vertices in the middle of placer centroids.

        Normally, linear vertices are in the middle of the
        non-linear vertices they are neighbors to. However,
        they should really be placed in the middle of the functional
        groups sitting on the vertices. This ensures that the bonds
        made are straight.

        For example, take ``N`` to be a non-linear vertex, ``L``
        to be a linear vertex and ``X`` to be a functional group
        placed on a vertex. If, vertex positions are not updated,
        a linear vertex will be in the middle of the two non-linear
        vertices::

            X               X
            |               |
            N    X--L--X    N

        Therefore, the bonds made will be bent::

            X           X
            | \       / |
            N  X--L--X  N

        However, this method updates the position of ``L`` to
        be in the middle of the functional groups::

            X    X--L--X    X
            |               |
            N               N

        Now when bonds are made, they will be straight::

            X----X--L--X----X
            |               |
            N               N

        Note that all vertices, which can have their positions updated,
        will have them updated, even the non-linear ones. A vertex
        can have its position updated, if the positions of all of the
        neighboring functional groups (the ``X`` in the diagrams
        above), are known.

        Parameters:

            vertices:
                The vertices which were just used to place a set of
                `building_blocks`.

            edges:
                the :class:`.Edge` instances it is connected to.

            building_blocks:
                The building blocks which are just placed.

            results:
                For each vertex in `vertices`, the result of the
                placement.

        """

        for vertex, vertex_edges, building_block, result in zip(
            vertices,
            edges,
            building_blocks,
            results,
        ):
            building_block = building_block.with_position_matrix(
                position_matrix=result.position_matrix,
            )
            edge_functional_groups = dict(zip(
                result.functional_group_edges.values(),
                result.functional_group_edges.keys(),
            ))
            for neighbor_id, edge_id in self._get_neighbors(
                vertex=vertex,
                vertex_edges=vertex_edges,
            ):
                fg_id = edge_functional_groups[edge_id]
                functional_group, = (
                    building_block.get_functional_groups(fg_id)
                )
                self._neighbor_positions[neighbor_id] = (
                    self._neighbor_positions.get(neighbor_id, [])
                )
                self._neighbor_positions[neighbor_id].append(
                    building_block.get_centroid(
                        atom_ids=functional_group.get_placer_ids(),
                    )
                )

        self._graph_state = self._graph_state.with_vertices(
            vertices=self._get_new_vertices(),
        )

    def _get_neighbors(
        self,
        vertex: CageVertex,
        vertex_edges: abc.Collection[Edge],
    ) -> abc.Iterator[tuple[int, int]]:
        """
        Yield the neighbor vertices of `vertex`.

        Parameters:

            vertex:
                The vertex whose neighbors are desired.

            vertex_edges:
                The edges connected to `vertex`.

        Yields:

            The first element is the id of a neighbor vertex and
            the second element is the id of the edge through which
            it is connected.

        """

        for edge in vertex_edges:
            neighbor_id = (
                edge.get_vertex1_id()
                if vertex.get_id() != edge.get_vertex1_id()
                else edge.get_vertex2_id()
            )
            neighbor, = self._graph_state.get_vertices(neighbor_id)
            if neighbor.use_neighbor_placement():
                yield neighbor_id, edge.get_id()

    def _get_new_vertices(self) -> abc.Iterator[CageVertex]:
        """
        Yield the vertices once new positions have been added.

        Yields:

            A vertex of the topology graph.

        """

        for vertex_id, vertex in enumerate(
            self._graph_state.get_vertices(),
        ):
            if (
                len(self._neighbor_positions.get(vertex_id, []))
                == self._vertex_degrees[vertex_id]
            ):

                neighbor_positions = self._neighbor_positions[
                    vertex_id
                ]
                position: np.ndarray = (
                    # I'm sure this is a ndarray.
                    sum(neighbor_positions)  # type: ignore
                    / len(neighbor_positions)
                )
                yield vertex.with_position(position)
            else:
                yield vertex

    def _clone(self: _T) -> _T:
        clone = super()._clone()
        clone._neighbor_positions = {
            key: list(value)
            for key, value in self._neighbor_positions.items()
        }
        clone._num_placement_stages_done = (
            self._num_placement_stages_done
        )
        clone._num_placement_stages = self._num_placement_stages
        clone._vertex_degrees = dict(self._vertex_degrees)
        return clone

    def clone(self) -> CageConstructionState:
        return self._clone()
