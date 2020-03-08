import numpy as np

from ..topology_graph import Vertex


class _HostVertex(Vertex):
    """
    Places the host in a :class:`.Complex`.

    """

    def place_building_block(self, building_block, edges):
        return building_block.with_centroid(
            position=self._position,
            atom_ids=building_block.get_placer_ids(),
        ).get_position_matrix()

    def map_functional_groups_to_edges(self, building_block, edges):
        return {}


class _GuestVertex(Vertex):
    """
    Places the guest in a :class:`.Complex`.

    """

    def __init__(self, id, position, start, target):
        self._start = np.array(start, dtype=np.float64)
        self._target = np.array(target, dtype=np.float64)
        super().__init__(id, position)

    def clone(self):
        clone = super().clone()
        clone._start = np.array(self._start)
        clone._target = np.array(self._target)
        return clone

    def place_building_block(self, building_block, edges):
        building_block = building_block.with_centroid(
            posiion=self._position,
            atom_ids=building_block.get_placer_ids(),
        )
        return building_block.with_rotation_between_vectors(
            start=self._start,
            target=self._target,
            origin=self._position,
        ).get_position_matrix()

    def map_functional_groups_to_edges(self, building_block, edges):
        return {}
