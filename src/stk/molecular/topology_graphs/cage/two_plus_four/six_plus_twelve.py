"""
Six Plus Twelve
===============

"""

from ..cage import Cage
from ..vertices import _LinearCageVertex, _NonLinearCageVertex
from ...topology_graph import Edge


class SixPlusTwelve(Cage):
    """
    Represents a cage topology graph.

    See :class:`.Cage` for more details and examples.

    """

    _vertex_prototypes = (
        _NonLinearCageVertex(0, [-1, -1, 0]),
        _NonLinearCageVertex(1, [-1, 1, 0]),
        _NonLinearCageVertex(2, [1, -1, 0]),
        _NonLinearCageVertex(3, [1, 1, 0]),
        _NonLinearCageVertex(4, [0, 0, 1]),
        _NonLinearCageVertex(5, [0, 0, -1]),
    )

    _vertex_prototypes = (
        *_vertex_prototypes,
        _LinearCageVertex.init_at_center(
            id=6,
            vertices=(_vertex_prototypes[0], _vertex_prototypes[1]),
        ),
        _LinearCageVertex.init_at_center(
            id=7,
            vertices=(_vertex_prototypes[1], _vertex_prototypes[3]),
        ),
        _LinearCageVertex.init_at_center(
            id=8,
            vertices=(_vertex_prototypes[3], _vertex_prototypes[2]),
        ),

        _LinearCageVertex.init_at_center(
            id=9,
            vertices=(_vertex_prototypes[0], _vertex_prototypes[2]),
        ),
        _LinearCageVertex.init_at_center(
            id=10,
            vertices=(_vertex_prototypes[4], _vertex_prototypes[0]),
        ),
        _LinearCageVertex.init_at_center(
            id=11,
            vertices=(_vertex_prototypes[4], _vertex_prototypes[1]),
        ),

        _LinearCageVertex.init_at_center(
            id=12,
            vertices=(_vertex_prototypes[4], _vertex_prototypes[2]),
        ),
        _LinearCageVertex.init_at_center(
            id=13,
            vertices=(_vertex_prototypes[4], _vertex_prototypes[3]),
        ),
        _LinearCageVertex.init_at_center(
            id=14,
            vertices=(_vertex_prototypes[5], _vertex_prototypes[0]),
        ),

        _LinearCageVertex.init_at_center(
            id=15,
            vertices=(_vertex_prototypes[5], _vertex_prototypes[1]),
        ),
        _LinearCageVertex.init_at_center(
            id=16,
            vertices=(_vertex_prototypes[5], _vertex_prototypes[2]),
        ),
        _LinearCageVertex.init_at_center(
            id=17,
            vertices=(_vertex_prototypes[5], _vertex_prototypes[3]),
        )
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[6], _vertex_prototypes[0]),
        Edge(1, _vertex_prototypes[6], _vertex_prototypes[1]),
        Edge(2, _vertex_prototypes[7], _vertex_prototypes[1]),
        Edge(3, _vertex_prototypes[7], _vertex_prototypes[3]),
        Edge(4, _vertex_prototypes[8], _vertex_prototypes[3]),
        Edge(5, _vertex_prototypes[8], _vertex_prototypes[2]),

        Edge(6, _vertex_prototypes[9], _vertex_prototypes[0]),
        Edge(7, _vertex_prototypes[9], _vertex_prototypes[2]),
        Edge(8, _vertex_prototypes[10], _vertex_prototypes[4]),
        Edge(9, _vertex_prototypes[10], _vertex_prototypes[0]),
        Edge(10, _vertex_prototypes[11], _vertex_prototypes[4]),
        Edge(11, _vertex_prototypes[11], _vertex_prototypes[1]),

        Edge(12, _vertex_prototypes[12], _vertex_prototypes[4]),
        Edge(13, _vertex_prototypes[12], _vertex_prototypes[2]),
        Edge(14, _vertex_prototypes[13], _vertex_prototypes[4]),
        Edge(15, _vertex_prototypes[13], _vertex_prototypes[3]),
        Edge(16, _vertex_prototypes[14], _vertex_prototypes[5]),
        Edge(17, _vertex_prototypes[14], _vertex_prototypes[0]),

        Edge(18, _vertex_prototypes[15], _vertex_prototypes[5]),
        Edge(19, _vertex_prototypes[15], _vertex_prototypes[1]),
        Edge(20, _vertex_prototypes[16], _vertex_prototypes[5]),
        Edge(21, _vertex_prototypes[16], _vertex_prototypes[2]),
        Edge(22, _vertex_prototypes[17], _vertex_prototypes[5]),
        Edge(23, _vertex_prototypes[17], _vertex_prototypes[3]),
    )

    _num_windows = 8
    _num_window_types = 1
