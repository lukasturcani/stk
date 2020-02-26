from typing import NamedTuple


class _PlacementResult(NamedTuple):
    """
    The result of a building block placement.

    Attributes
    ----------
    position_matrix : :class:`numpy.ndarray`
        The position matrix of the building block after it has been
        placed on the vertex.

    functional_group_edges : :class:`dict`
        Maps the id of a functional group to the id of the edge it is
        assigned to.

    """

    position_matrix: object
    functional_group_edges: object


class _Placement:
    def __init__(self, vertex, edges, building_block):
        self._vertex = vertex
        self._edges = edges
        self._building_block = building_block

    def get_result(self):
        position_matrix = self._vertex.place_building_block(
            building_block=self._building_block,
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
