"""
Macrocycle Vertices
===================

"""

from __future__ import annotations

import typing

import numpy as np
from scipy.spatial.distance import euclidean

from stk._internal.building_block import BuildingBlock
from stk._internal.topology_graphs.edge import Edge
from stk._internal.topology_graphs.vertex import Vertex


class CycleVertex(Vertex):
    """
    Represents a vertex in a macrocycle.

    """

    def __init__(
        self,
        id: int,
        position: typing.Union[np.ndarray, tuple[float, float, float]],
        flip: bool,
        angle: float,
    ) -> None:
        """
        Initialize a :class:`.CycleVertex` instance.

        Parameters:

            id:
                The id of the vertex.

            position:
                The position of the vertex.

            flip:
                If ``True``, the orientation of building blocks placed
                by the vertex will be flipped.

            angle:
                The position of the vertex along the cycle, specified
                by the `angle`.

        """

        super().__init__(id, position)
        self._flip = flip
        self._angle = angle

    def clone(self) -> CycleVertex:
        clone = self._clone()
        clone._flip = self._flip
        clone._angle = self._angle
        return clone

    def get_flip(self) -> bool:
        """
        Return ``True`` if the vertex flips building blocks it places.

        Returns:

            ``True`` if the vertex flips building blocks it places.

        """

        return self._flip

    def place_building_block(
        self,
        building_block: BuildingBlock,
        edges: tuple[Edge, ...],
    ) -> np.ndarray:
        assert building_block.get_num_functional_groups() == 2, (
            f"{building_block} needs to have exactly 2 functional "
            "groups but has "
            f"{building_block.get_num_functional_groups()}."
        )
        building_block = building_block.with_centroid(
            position=self._position,
            atom_ids=building_block.get_placer_ids(),
        )
        fg0, fg1 = building_block.get_functional_groups()
        fg0_position = building_block.get_centroid(
            atom_ids=fg0.get_placer_ids(),
        )
        fg1_position = building_block.get_centroid(
            atom_ids=fg1.get_placer_ids(),
        )
        return (
            building_block.with_rotation_between_vectors(
                start=fg1_position - fg0_position,
                target=np.array([-1 if self._flip else 1, 0, 0]),
                origin=self._position,
            )
            .with_rotation_about_axis(
                angle=self._angle - (np.pi / 2),
                axis=np.array([0, 0, 1]),
                origin=self._position,
            )
            .get_position_matrix()
        )

    def map_functional_groups_to_edges(
        self,
        building_block: BuildingBlock,
        edges: tuple[Edge, ...],
    ) -> dict[int, int]:
        fg0_position = building_block.get_centroid(
            atom_ids=next(
                building_block.get_functional_groups()
            ).get_placer_ids(),
        )

        def fg0_distance(edge: Edge) -> float:
            return euclidean(edge.get_position(), fg0_position)

        edge0 = min(edges, key=fg0_distance)
        return {
            0: edge0.get_id(),
            1: (edges[1].get_id() if edge0 is edges[0] else edges[0].get_id()),
        }

    def __str__(self) -> str:
        return (
            f"Vertex(id={self._id}, "
            f"position={tuple(self._position.tolist())}, "
            f"flip={self._flip}, "
            f"angle={self._angle})"
        )
