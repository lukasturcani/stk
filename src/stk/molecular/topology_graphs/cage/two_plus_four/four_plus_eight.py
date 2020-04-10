"""
Four Plus Eight
===============

"""

from ..cage import Cage
from ..vertices import _LinearCageVertex, _NonLinearCageVertex
from ...topology_graph import Edge


class FourPlusEight(Cage):
    """
    Represents a cage topology graph.

    See :class:`.Cage` for more details and examples.

    """

    _vertex_prototypes = (
        _NonLinearCageVertex(0, [-1, -1, 0], False),
        _NonLinearCageVertex(1, [-1, 1, 0], False),

        _NonLinearCageVertex(2, [1, -1, 0], False),
        _NonLinearCageVertex(3, [1, 1, 0], False),

        _LinearCageVertex(4, [-2, 0, 1], False),
        _LinearCageVertex(5, [-2, 0, -1], False),

        _LinearCageVertex(6, [0, 2, 1], False),
        _LinearCageVertex(7, [0, 2, -1], False),

        _LinearCageVertex(8, [0, -2, 1], False),
        _LinearCageVertex(9, [0, -2, -1], False),

        _LinearCageVertex(10, [2, 0, 1], False),
        _LinearCageVertex(11, [2, 0, -1], False),

    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[4], _vertex_prototypes[0]),
        Edge(1, _vertex_prototypes[4], _vertex_prototypes[1]),

        Edge(2, _vertex_prototypes[5], _vertex_prototypes[0]),
        Edge(3, _vertex_prototypes[5], _vertex_prototypes[1]),

        Edge(4, _vertex_prototypes[6], _vertex_prototypes[1]),
        Edge(5, _vertex_prototypes[6], _vertex_prototypes[3]),

        Edge(6, _vertex_prototypes[7], _vertex_prototypes[1]),
        Edge(7, _vertex_prototypes[7], _vertex_prototypes[3]),

        Edge(8, _vertex_prototypes[8], _vertex_prototypes[0]),
        Edge(9, _vertex_prototypes[8], _vertex_prototypes[2]),

        Edge(10, _vertex_prototypes[9], _vertex_prototypes[0]),
        Edge(11, _vertex_prototypes[9], _vertex_prototypes[2]),

        Edge(12, _vertex_prototypes[10], _vertex_prototypes[2]),
        Edge(13, _vertex_prototypes[10], _vertex_prototypes[3]),

        Edge(14, _vertex_prototypes[11], _vertex_prototypes[2]),
        Edge(15, _vertex_prototypes[11], _vertex_prototypes[3]),

    )

    _num_windows = 6
    _num_window_types = 2
