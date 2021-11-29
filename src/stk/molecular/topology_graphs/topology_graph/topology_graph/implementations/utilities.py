"""
Topology Graph Implementation Utilities
=======================================

"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np


class _PlacementResult(NamedTuple):
    """
    The result of a building block placement.

    Attributes:

        position_matrix:
            The position matrix of the building block after it has been
            placed on the vertex.

        functional_group_edges:
            Maps the id of a functional group to the id of the edge it
            is assigned to.

    """

    position_matrix: np.ndarray
    functional_group_edges: dict[int, int]


class _Placement:
    """
    Represents placement of a building block on a vertex.

    It represents a computation which carries out the placement
    of the building block on the vertex and the mapping of its
    functional group to edges.

    """

    def __init__(self, vertex, edges, building_block):
        """
        Initialize a :class:`._Placement`.

        Parameters
        ----------
        vertex : :class:`.Vertex`
            The vertex which does the placement.

        edges : :class:`tuple` of :class:`.Edge`
            The edges connected to `vertex`.

        building_block : :class:`.BuildingBlock`
            The building block to be placed on `vertex`.

        """

        self._vertex = vertex
        self._edges = edges
        self._building_block = building_block

    def get_result(self):
        """
        Get the result of the placement.

        Returns
        -------
        :class:`_PlacementResult`
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
        return _PlacementResult(
            position_matrix=position_matrix,
            functional_group_edges=functional_group_edges,
        )
