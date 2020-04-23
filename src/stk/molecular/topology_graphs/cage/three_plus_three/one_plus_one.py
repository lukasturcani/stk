"""
One Plus One
============

"""

import numpy as np

from stk.utilities import get_acute_vector
from ..vertices import _NonLinearCageVertex
from ..cage import Cage
from ...topology_graph import Edge


class _OnePlusOneVertex(_NonLinearCageVertex):
    def __init__(
        self,
        id,
        position,
        edge_normal,
        use_neighbor_placement=True,
        aligner_edge=0,
    ):
        super().__init__(
            id=id,
            position=position,
            use_neighbor_placement=use_neighbor_placement,
            aligner_edge=aligner_edge,
        )
        self._edge_normal = np.array(edge_normal, dtype=np.float64)

    def clone(self):
        clone = super().clone()
        clone._edge_normal = np.array(self._edge_normal)
        return clone

    def place_building_block(self, building_block, edges):
        building_block = building_block.with_centroid(
            position=self._position,
            atom_ids=building_block.get_placer_ids(),
        )
        edge_centroid = (
            sum(edge.get_position() for edge in edges) / len(edges)
        )
        core_centroid = building_block.get_centroid(
            atom_ids=building_block.get_core_atom_ids(),
        )
        placer_centroid = building_block.get_centroid(
            atom_ids=building_block.get_placer_ids(),
        )
        building_block = building_block.with_rotation_between_vectors(
            start=get_acute_vector(
                reference=core_centroid - placer_centroid,
                vector=building_block.get_plane_normal(
                    atom_ids=building_block.get_placer_ids(),
                ),
            ),
            target=self._edge_normal,
            origin=self._position,
        )
        fg_bonder_centroid = building_block.get_centroid(
            atom_ids=next(
                building_block.get_functional_groups()
            ).get_placer_ids(),
        )
        start = fg_bonder_centroid - self._position
        edge_coord = edges[self._aligner_edge].get_position()
        building_block = (
            building_block.with_rotation_to_minimize_angle(
                start=start,
                target=edge_coord - edge_centroid,
                axis=self._edge_normal,
                origin=self._position,
            )
        )
        return building_block.get_position_matrix()


class OnePlusOne(Cage):
    """
    Represents a capsule cage topology graph.

    See :class:`.Cage` for more details and examples.

    """

    _x = 1
    _vertex_prototypes = (
        _OnePlusOneVertex(0, [_x, 0., 0.], [1, 0, 0], False),
        _OnePlusOneVertex(1, [-_x, 0., 0.], [-1, 0, 0], False),

    )
    _edge_prototypes = (
        Edge(
            id=0,
            vertex1=_vertex_prototypes[0],
            vertex2=_vertex_prototypes[1],
            position=np.array([0., 1., 0.]),
        ),
        Edge(
            id=1,
            vertex1=_vertex_prototypes[0],
            vertex2=_vertex_prototypes[1],
            position=np.array([0., -1., 1.]),
        ),
        Edge(
            id=2,
            vertex1=_vertex_prototypes[0],
            vertex2=_vertex_prototypes[1],
            position=np.array([0., -1., -1.]),
        ),
    )

    _num_windows = 3
    _num_window_types = 1
