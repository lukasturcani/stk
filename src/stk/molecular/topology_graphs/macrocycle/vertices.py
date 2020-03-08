import numpy as np
from functools import partial
from scipy.spatial.distance import euclidean

from ..topology_graph import Vertex


class _CycleVertex(Vertex):
    """
    Represents a vertex in a macrocycle.

    """

    def __init__(self, id, position, flip, angle):
        super().__init__(id, position)
        self._flip = flip
        self._angle = angle

    def clone(self):
        clone = super().clone()
        clone._flip = self._flip
        clone._angle = self._angle
        return clone

    def get_flip(self):
        return self._flip

    def place_building_block(self, building_block, edges):
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
        return building_block.with_rotation_between_vectors(
            start=fg1_position - fg0_position,
            target=[-1 if self._flip else 1, 0, 0],
            origin=self._position,
        ).with_rotation_about_axis(
            angle=self._angle-(np.pi/2),
            axis=np.array([0, 0, 1]),
            origin=self._position,
        ).get_position_matrix()

    def map_functional_groups_to_edges(self, building_block, edges):
        fg0_position = building_block.get_centroid(
            atom_ids=next(
                building_block.get_functional_groups()
            ).get_placer_ids(),
        )
        fg0_distance = partial(euclidean, fg0_position)
        edge0 = min(
            edges,
            key=lambda edge: fg0_distance(edge.get_position()),
        )
        return {
            0: edge0.get_id(),
            1: (
                edges[1].get_id()
                if edge0 is edges[0]
                else edges[0].get_id()
            ),
        }

    def __str__(self):
        return (
            f'Vertex(id={self.id}, '
            f'position={self._position.tolist()}, '
            f'flip={self._flip}, '
            f'angle={self._angle})'
        )
