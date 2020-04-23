"""
Two Plus Two
============

"""

import numpy as np

from ..cage import Cage
from ..vertices import _NonLinearCageVertex
from ...topology_graph import Edge


class TwoPlusTwo(Cage):
    """
    Represents a tetrahedron cage topology graph.

    See :class:`.Cage` for more details and examples.

    """

    _x = 1
    _vertex_prototypes = (
        _NonLinearCageVertex(0, [_x, 0, -_x/np.sqrt(2)], False),
        _NonLinearCageVertex(1, [-_x, 0, -_x/np.sqrt(2)], False),
        _NonLinearCageVertex(2, [0, _x, _x/np.sqrt(2)], False),
        _NonLinearCageVertex(3, [0, -_x, _x/np.sqrt(2)], False)
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[1]),
        Edge(1, _vertex_prototypes[0], _vertex_prototypes[2]),
        Edge(2, _vertex_prototypes[0], _vertex_prototypes[3]),

        Edge(3, _vertex_prototypes[1], _vertex_prototypes[2]),
        Edge(4, _vertex_prototypes[1], _vertex_prototypes[3]),

        Edge(5, _vertex_prototypes[2], _vertex_prototypes[3])
    )

    _num_windows = 4
    _num_window_types = 1
