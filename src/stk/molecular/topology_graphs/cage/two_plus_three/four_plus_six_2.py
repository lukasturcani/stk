"""
Four Plus Six 2
===============

"""

from ..cage import Cage
from ..vertices import _LinearCageVertex, _NonLinearCageVertex
from ...topology_graph import Edge


class FourPlusSix2(Cage):
    """
    Represents a cage topology graph.

    See :class:`.Cage` for more details and examples.

    """

    _vertex_prototypes = (
        _NonLinearCageVertex(0, [1, 0, 1]),
        _NonLinearCageVertex(1, [-1, 0, 1]),
        _NonLinearCageVertex(2, [1, 0, -1]),
        _NonLinearCageVertex(3, [-1, 0, -1]),

        _LinearCageVertex(4, [0, -1, 1], False),
        _LinearCageVertex(5, [0, 1, 1], False),
        _LinearCageVertex(6, [0, -1, -1], False),
        _LinearCageVertex(7, [0, 1, -1], False),
    )

    _vertex_prototypes = (
        *_vertex_prototypes,
        _LinearCageVertex.init_at_center(
            id=8,
            vertices=(_vertex_prototypes[0], _vertex_prototypes[2]),
        ),
        _LinearCageVertex.init_at_center(
            id=9,
            vertices=(_vertex_prototypes[1], _vertex_prototypes[3]),
        )
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[4], _vertex_prototypes[0]),
        Edge(1, _vertex_prototypes[4], _vertex_prototypes[1]),
        Edge(2, _vertex_prototypes[5], _vertex_prototypes[0]),
        Edge(3, _vertex_prototypes[5], _vertex_prototypes[1]),
        Edge(4, _vertex_prototypes[6], _vertex_prototypes[2]),
        Edge(5, _vertex_prototypes[6], _vertex_prototypes[3]),
        Edge(6, _vertex_prototypes[7], _vertex_prototypes[2]),
        Edge(7, _vertex_prototypes[7], _vertex_prototypes[3]),
        Edge(8, _vertex_prototypes[8], _vertex_prototypes[0]),
        Edge(9, _vertex_prototypes[8], _vertex_prototypes[2]),
        Edge(10, _vertex_prototypes[9], _vertex_prototypes[1]),
        Edge(11, _vertex_prototypes[9], _vertex_prototypes[3]),
    )

    _num_windows = 4
    _num_window_types = 2
