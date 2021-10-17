"""
Topology Graph Implementation Utilities
=======================================

"""

from __future__ import annotations

from collections import abc

from .....building_block import BuildingBlock
from ...placement_result import PlacementResult
from ...vertex import Vertex
from ...edge import Edge


__all__ = (
    'Placement',
)


class Placement:
    """
    Represents placement of a building block on a vertex.

    It represents a computation which carries out the placement
    of the building block on the vertex and the mapping of its
    functional group to edges.

    """

    def __init__(
        self,
        vertex: Vertex,
        edges: abc.Collection[Edge],
        building_block: BuildingBlock,
    ) -> None:
        """
        Initialize a :class:`.Placement`.

        Parameters:

            vertex:
                The vertex which does the placement.

            edges:
                The edges connected to `vertex`.

            building_block:
                The building block to be placed on `vertex`.

        """

        self._vertex = vertex
        self._edges = edges
        self._building_block = building_block

    def get_result(self) -> PlacementResult:
        """
        Get the result of the placement.

        Returns:

            The result of the placement.

        """

        position_matrix = self._vertex.place_building_block(
            building_block=self._building_block,
            edges=self._edges,
        )
        position_matrix.setflags(write=False)
        building_block = self._building_block.with_position_matrix(
            position_matrix=position_matrix,
        )
        functional_group_edges = (
            self._vertex.map_functional_groups_to_edges(
                building_block=building_block,
                edges=self._edges,
            )
        )
        return PlacementResult(
            position_matrix=position_matrix,
            functional_group_edges=functional_group_edges,
        )
