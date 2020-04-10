"""
Ten Plus Twenty
===============

"""

from ..cage import Cage
from ..vertices import _LinearCageVertex, _NonLinearCageVertex
from ...topology_graph import Edge


class TenPlusTwenty(Cage):
    """
    Represents a cage topology graph.

    See :class:`.Cage` for more details and examples.

    """

    _x = 1.5
    _z = _x / 2
    _vertex_prototypes = (
        _NonLinearCageVertex(0, [-_x, _x, -_z]),
        _NonLinearCageVertex(1, [-_x, -_x, -_z]),
        _NonLinearCageVertex(2, [_x, _x, -_z]),
        _NonLinearCageVertex(3, [_x, -_x, -_z]),

        _NonLinearCageVertex(4, [-_x, _x, _z]),
        _NonLinearCageVertex(5, [-_x, -_x, _z]),
        _NonLinearCageVertex(6, [_x, _x, _z]),
        _NonLinearCageVertex(7, [_x, -_x, _z]),

        _NonLinearCageVertex(8, [0, 0, _x]),
        _NonLinearCageVertex(9, [0, 0, -_x]),
    )

    _vertex_prototypes = (
        *_vertex_prototypes,

        _LinearCageVertex.init_at_center(
            id=10,
            vertices=(_vertex_prototypes[0], _vertex_prototypes[2]),
        ),
        _LinearCageVertex.init_at_center(
            id=11,
            vertices=(_vertex_prototypes[0], _vertex_prototypes[1]),
        ),
        _LinearCageVertex.init_at_center(
            id=12,
            vertices=(_vertex_prototypes[1], _vertex_prototypes[3]),
        ),
        _LinearCageVertex.init_at_center(
            id=13,
            vertices=(_vertex_prototypes[2], _vertex_prototypes[3]),
        ),

        _LinearCageVertex.init_at_center(
            id=14,
            vertices=(_vertex_prototypes[4], _vertex_prototypes[6]),
        ),
        _LinearCageVertex.init_at_center(
            id=15,
            vertices=(_vertex_prototypes[4], _vertex_prototypes[5]),
        ),
        _LinearCageVertex.init_at_center(
            id=16,
            vertices=(_vertex_prototypes[5], _vertex_prototypes[7]),
        ),
        _LinearCageVertex.init_at_center(
            id=17,
            vertices=(_vertex_prototypes[6], _vertex_prototypes[7]),
        ),

        _LinearCageVertex.init_at_center(
            id=18,
            vertices=(_vertex_prototypes[0], _vertex_prototypes[4]),
        ),
        _LinearCageVertex.init_at_center(
            id=19,
            vertices=(_vertex_prototypes[1], _vertex_prototypes[5]),
        ),
        _LinearCageVertex.init_at_center(
            id=20,
            vertices=(_vertex_prototypes[2], _vertex_prototypes[6]),
        ),
        _LinearCageVertex.init_at_center(
            id=21,
            vertices=(_vertex_prototypes[3], _vertex_prototypes[7]),
        ),

        _LinearCageVertex.init_at_center(
            id=22,
            vertices=(_vertex_prototypes[8], _vertex_prototypes[4]),
        ),
        _LinearCageVertex.init_at_center(
            id=23,
            vertices=(_vertex_prototypes[8], _vertex_prototypes[5]),
        ),
        _LinearCageVertex.init_at_center(
            id=24,
            vertices=(_vertex_prototypes[8], _vertex_prototypes[6]),
        ),
        _LinearCageVertex.init_at_center(
            id=25,
            vertices=(_vertex_prototypes[8], _vertex_prototypes[7]),
        ),

        _LinearCageVertex.init_at_center(
            id=26,
            vertices=(_vertex_prototypes[9], _vertex_prototypes[0]),
        ),
        _LinearCageVertex.init_at_center(
            id=27,
            vertices=(_vertex_prototypes[9], _vertex_prototypes[1]),
        ),
        _LinearCageVertex.init_at_center(
            id=28,
            vertices=(_vertex_prototypes[9], _vertex_prototypes[2]),
        ),
        _LinearCageVertex.init_at_center(
            id=29,
            vertices=(_vertex_prototypes[9], _vertex_prototypes[3]),
        )
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[10], _vertex_prototypes[0]),
        Edge(1, _vertex_prototypes[10], _vertex_prototypes[2]),
        Edge(2, _vertex_prototypes[11], _vertex_prototypes[0]),
        Edge(3, _vertex_prototypes[11], _vertex_prototypes[1]),

        Edge(4, _vertex_prototypes[12], _vertex_prototypes[1]),
        Edge(5, _vertex_prototypes[12], _vertex_prototypes[3]),
        Edge(6, _vertex_prototypes[13], _vertex_prototypes[2]),
        Edge(7, _vertex_prototypes[13], _vertex_prototypes[3]),

        Edge(8, _vertex_prototypes[14], _vertex_prototypes[4]),
        Edge(9, _vertex_prototypes[14], _vertex_prototypes[6]),
        Edge(10, _vertex_prototypes[15], _vertex_prototypes[4]),
        Edge(11, _vertex_prototypes[15], _vertex_prototypes[5]),

        Edge(12, _vertex_prototypes[16], _vertex_prototypes[5]),
        Edge(13, _vertex_prototypes[16], _vertex_prototypes[7]),
        Edge(14, _vertex_prototypes[17], _vertex_prototypes[6]),
        Edge(15, _vertex_prototypes[17], _vertex_prototypes[7]),

        Edge(16, _vertex_prototypes[18], _vertex_prototypes[0]),
        Edge(17, _vertex_prototypes[18], _vertex_prototypes[4]),
        Edge(18, _vertex_prototypes[19], _vertex_prototypes[1]),
        Edge(19, _vertex_prototypes[19], _vertex_prototypes[5]),

        Edge(20, _vertex_prototypes[20], _vertex_prototypes[2]),
        Edge(21, _vertex_prototypes[20], _vertex_prototypes[6]),
        Edge(22, _vertex_prototypes[21], _vertex_prototypes[3]),
        Edge(23, _vertex_prototypes[21], _vertex_prototypes[7]),

        Edge(24, _vertex_prototypes[22], _vertex_prototypes[8]),
        Edge(25, _vertex_prototypes[22], _vertex_prototypes[4]),
        Edge(26, _vertex_prototypes[23], _vertex_prototypes[8]),
        Edge(27, _vertex_prototypes[23], _vertex_prototypes[5]),

        Edge(28, _vertex_prototypes[24], _vertex_prototypes[8]),
        Edge(29, _vertex_prototypes[24], _vertex_prototypes[6]),
        Edge(30, _vertex_prototypes[25], _vertex_prototypes[8]),
        Edge(31, _vertex_prototypes[25], _vertex_prototypes[7]),

        Edge(32, _vertex_prototypes[26], _vertex_prototypes[9]),
        Edge(33, _vertex_prototypes[26], _vertex_prototypes[0]),
        Edge(34, _vertex_prototypes[27], _vertex_prototypes[9]),
        Edge(35, _vertex_prototypes[27], _vertex_prototypes[1]),

        Edge(36, _vertex_prototypes[28], _vertex_prototypes[9]),
        Edge(37, _vertex_prototypes[28], _vertex_prototypes[2]),
        Edge(38, _vertex_prototypes[29], _vertex_prototypes[9]),
        Edge(39, _vertex_prototypes[29], _vertex_prototypes[3]),
    )

    _num_windows = 12
    _num_window_types = 2
