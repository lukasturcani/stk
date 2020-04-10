"""
Twelve Plus Thirty
==================

"""

from scipy.constants import golden

from ..cage import Cage
from ..vertices import _LinearCageVertex, _NonLinearCageVertex
from ...topology_graph import Edge


class TwelvePlusThirty(Cage):
    """
    Represents a icosahedron cage topology graph.

    See :class:`.Cage` for more details and examples.

    """
    # Vertices of a regular origin-centred icosahedron
    # Source: http://eusebeia.dyndns.org/4d/icosahedron
    _vertex_prototypes = (
        _NonLinearCageVertex(0, [0, 1, golden]),
        _NonLinearCageVertex(1, [0, -1, golden]),
        _NonLinearCageVertex(2, [0, 1, -golden]),
        _NonLinearCageVertex(3, [0, -1, -golden]),
        _NonLinearCageVertex(4, [1, golden, 0]),
        _NonLinearCageVertex(5, [-1, golden, 0]),
        _NonLinearCageVertex(6, [1, -golden, 0]),
        _NonLinearCageVertex(7, [-1, -golden, 0]),
        _NonLinearCageVertex(8, [golden, 0, 1]),
        _NonLinearCageVertex(9, [-golden, 0, 1]),
        _NonLinearCageVertex(10, [golden, 0, -1]),
        _NonLinearCageVertex(11, [-golden, 0, -1])
    )

    _vertex_prototypes = (
        *_vertex_prototypes,
        _LinearCageVertex.init_at_center(
            id=12,
            vertices=(_vertex_prototypes[0], _vertex_prototypes[1]),
        ),
        _LinearCageVertex.init_at_center(
            id=13,
            vertices=(_vertex_prototypes[0], _vertex_prototypes[9]),
        ),
        _LinearCageVertex.init_at_center(
            id=14,
            vertices=(_vertex_prototypes[0], _vertex_prototypes[5]),
        ),
        _LinearCageVertex.init_at_center(
            id=15,
            vertices=(_vertex_prototypes[0], _vertex_prototypes[4]),
        ),
        _LinearCageVertex.init_at_center(
            id=16,
            vertices=(_vertex_prototypes[0], _vertex_prototypes[8]),
        ),
        _LinearCageVertex.init_at_center(
            id=17,
            vertices=(_vertex_prototypes[8], _vertex_prototypes[1]),),
        _LinearCageVertex.init_at_center(
            id=18,
            vertices=(_vertex_prototypes[1], _vertex_prototypes[9]),
        ),
        _LinearCageVertex.init_at_center(
            id=19,
            vertices=(_vertex_prototypes[9], _vertex_prototypes[5]),
        ),
        _LinearCageVertex.init_at_center(
            id=20,
            vertices=(_vertex_prototypes[5], _vertex_prototypes[4]),
        ),
        _LinearCageVertex.init_at_center(
            id=21,
            vertices=(_vertex_prototypes[4], _vertex_prototypes[8]),
        ),
        _LinearCageVertex.init_at_center(
            id=22,
            vertices=(_vertex_prototypes[5], _vertex_prototypes[2]),
        ),
        _LinearCageVertex.init_at_center(
            id=23,
            vertices=(_vertex_prototypes[5], _vertex_prototypes[11]),
        ),
        _LinearCageVertex.init_at_center(
            id=24,
            vertices=(_vertex_prototypes[9], _vertex_prototypes[11]),
        ),
        _LinearCageVertex.init_at_center(
            id=25,
            vertices=(_vertex_prototypes[9], _vertex_prototypes[7]),
        ),
        _LinearCageVertex.init_at_center(
            id=26,
            vertices=(_vertex_prototypes[1], _vertex_prototypes[7]),
        ),
        _LinearCageVertex.init_at_center(
            id=27,
            vertices=(_vertex_prototypes[1], _vertex_prototypes[6]),
        ),
        _LinearCageVertex.init_at_center(
            id=28,
            vertices=(_vertex_prototypes[8], _vertex_prototypes[6]),
        ),
        _LinearCageVertex.init_at_center(
            id=29,
            vertices=(_vertex_prototypes[8], _vertex_prototypes[10]),
        ),
        _LinearCageVertex.init_at_center(
            id=30,
            vertices=(_vertex_prototypes[4], _vertex_prototypes[10]),
        ),
        _LinearCageVertex.init_at_center(
            id=31,
            vertices=(_vertex_prototypes[4], _vertex_prototypes[2]),
        ),
        _LinearCageVertex.init_at_center(
            id=32,
            vertices=(_vertex_prototypes[2], _vertex_prototypes[11]),
        ),
        _LinearCageVertex.init_at_center(
            id=33,
            vertices=(_vertex_prototypes[11], _vertex_prototypes[7]),
        ),
        _LinearCageVertex.init_at_center(
            id=34,
            vertices=(_vertex_prototypes[7], _vertex_prototypes[6]),
        ),
        _LinearCageVertex.init_at_center(
            id=35,
            vertices=(_vertex_prototypes[6], _vertex_prototypes[10]),
        ),
        _LinearCageVertex.init_at_center(
            id=36,
            vertices=(_vertex_prototypes[10], _vertex_prototypes[2]),
        ),
        _LinearCageVertex.init_at_center(
            id=37,
            vertices=(_vertex_prototypes[2], _vertex_prototypes[3]),
        ),
        _LinearCageVertex.init_at_center(
            id=38,
            vertices=(_vertex_prototypes[11], _vertex_prototypes[3]),
        ),
        _LinearCageVertex.init_at_center(
            id=39,
            vertices=(_vertex_prototypes[7], _vertex_prototypes[3]),
        ),
        _LinearCageVertex.init_at_center(
            id=40,
            vertices=(_vertex_prototypes[6], _vertex_prototypes[3]),
        ),
        _LinearCageVertex.init_at_center(
            id=41,
            vertices=(_vertex_prototypes[10], _vertex_prototypes[3]),
        )
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[12], _vertex_prototypes[0]),
        Edge(1, _vertex_prototypes[12], _vertex_prototypes[1]),
        Edge(2, _vertex_prototypes[13], _vertex_prototypes[0]),
        Edge(3, _vertex_prototypes[13], _vertex_prototypes[9]),

        Edge(4, _vertex_prototypes[14], _vertex_prototypes[0]),
        Edge(5, _vertex_prototypes[14], _vertex_prototypes[5]),
        Edge(6, _vertex_prototypes[15], _vertex_prototypes[0]),
        Edge(7, _vertex_prototypes[15], _vertex_prototypes[4]),

        Edge(8, _vertex_prototypes[16], _vertex_prototypes[0]),
        Edge(9, _vertex_prototypes[16], _vertex_prototypes[8]),
        Edge(10, _vertex_prototypes[17], _vertex_prototypes[8]),
        Edge(11, _vertex_prototypes[17], _vertex_prototypes[1]),

        Edge(12, _vertex_prototypes[18], _vertex_prototypes[1]),
        Edge(13, _vertex_prototypes[18], _vertex_prototypes[9]),
        Edge(14, _vertex_prototypes[19], _vertex_prototypes[9]),
        Edge(15, _vertex_prototypes[19], _vertex_prototypes[5]),

        Edge(16, _vertex_prototypes[20], _vertex_prototypes[5]),
        Edge(17, _vertex_prototypes[20], _vertex_prototypes[4]),
        Edge(18, _vertex_prototypes[21], _vertex_prototypes[4]),
        Edge(19, _vertex_prototypes[21], _vertex_prototypes[8]),

        Edge(20, _vertex_prototypes[22], _vertex_prototypes[5]),
        Edge(21, _vertex_prototypes[22], _vertex_prototypes[2]),
        Edge(22, _vertex_prototypes[23], _vertex_prototypes[5]),
        Edge(23, _vertex_prototypes[23], _vertex_prototypes[11]),

        Edge(24, _vertex_prototypes[24], _vertex_prototypes[9]),
        Edge(25, _vertex_prototypes[24], _vertex_prototypes[11]),
        Edge(26, _vertex_prototypes[25], _vertex_prototypes[9]),
        Edge(27, _vertex_prototypes[25], _vertex_prototypes[7]),

        Edge(28, _vertex_prototypes[26], _vertex_prototypes[1]),
        Edge(29, _vertex_prototypes[26], _vertex_prototypes[7]),
        Edge(30, _vertex_prototypes[27], _vertex_prototypes[1]),
        Edge(31, _vertex_prototypes[27], _vertex_prototypes[6]),

        Edge(32, _vertex_prototypes[28], _vertex_prototypes[8]),
        Edge(33, _vertex_prototypes[28], _vertex_prototypes[6]),
        Edge(34, _vertex_prototypes[29], _vertex_prototypes[8]),
        Edge(35, _vertex_prototypes[29], _vertex_prototypes[10]),

        Edge(36, _vertex_prototypes[30], _vertex_prototypes[4]),
        Edge(37, _vertex_prototypes[30], _vertex_prototypes[10]),
        Edge(38, _vertex_prototypes[31], _vertex_prototypes[4]),
        Edge(39, _vertex_prototypes[31], _vertex_prototypes[2]),

        Edge(40, _vertex_prototypes[32], _vertex_prototypes[2]),
        Edge(41, _vertex_prototypes[32], _vertex_prototypes[11]),
        Edge(42, _vertex_prototypes[33], _vertex_prototypes[11]),
        Edge(43, _vertex_prototypes[33], _vertex_prototypes[7]),

        Edge(44, _vertex_prototypes[34], _vertex_prototypes[7]),
        Edge(45, _vertex_prototypes[34], _vertex_prototypes[6]),
        Edge(46, _vertex_prototypes[35], _vertex_prototypes[6]),
        Edge(47, _vertex_prototypes[35], _vertex_prototypes[10]),

        Edge(48, _vertex_prototypes[36], _vertex_prototypes[10]),
        Edge(49, _vertex_prototypes[36], _vertex_prototypes[2]),
        Edge(50, _vertex_prototypes[37], _vertex_prototypes[2]),
        Edge(51, _vertex_prototypes[37], _vertex_prototypes[3]),

        Edge(52, _vertex_prototypes[38], _vertex_prototypes[11]),
        Edge(53, _vertex_prototypes[38], _vertex_prototypes[3]),
        Edge(54, _vertex_prototypes[39], _vertex_prototypes[7]),
        Edge(55, _vertex_prototypes[39], _vertex_prototypes[3]),

        Edge(56, _vertex_prototypes[40], _vertex_prototypes[6]),
        Edge(57, _vertex_prototypes[40], _vertex_prototypes[3]),
        Edge(58, _vertex_prototypes[41], _vertex_prototypes[10]),
        Edge(59, _vertex_prototypes[41], _vertex_prototypes[3]),
    )

    _num_windows = 20
    _num_window_types = 1
