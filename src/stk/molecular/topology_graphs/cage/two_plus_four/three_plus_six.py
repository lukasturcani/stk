"""
Three Plus Six
==============

"""

import numpy as np

from ..cage import Cage
from ..vertices import _LinearCageVertex, _NonLinearCageVertex
from ...topology_graph import Edge


class ThreePlusSix(Cage):
    """
    Represents a cage topology graph.

    See :class:`.Cage` for more details and examples.

    """

    _x = 1
    _vertex_prototypes = (
        _NonLinearCageVertex(0, [-2*_x, -_x*np.sqrt(3), 0], False),
        _NonLinearCageVertex(1, [2*_x, -_x*np.sqrt(3), 0], False),
        _NonLinearCageVertex(2, [0, _x*np.sqrt(3), 0], False),

        _LinearCageVertex(3, [0, -2*_x*np.sqrt(3), _x], False),
        _LinearCageVertex(4, [0, -2*_x*np.sqrt(3), -_x], False),

        _LinearCageVertex(5, [2*_x, 0, _x], False),
        _LinearCageVertex(6, [2*_x, 0, -_x], False),

        _LinearCageVertex(7, [-2*_x, 0, _x], False),
        _LinearCageVertex(8, [-2*_x, 0, -_x], False),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[3], _vertex_prototypes[0]),
        Edge(1, _vertex_prototypes[3], _vertex_prototypes[1]),

        Edge(2, _vertex_prototypes[4], _vertex_prototypes[0]),
        Edge(3, _vertex_prototypes[4], _vertex_prototypes[1]),

        Edge(4, _vertex_prototypes[5], _vertex_prototypes[1]),
        Edge(5, _vertex_prototypes[5], _vertex_prototypes[2]),

        Edge(6, _vertex_prototypes[6], _vertex_prototypes[1]),
        Edge(7, _vertex_prototypes[6], _vertex_prototypes[2]),

        Edge(8, _vertex_prototypes[7], _vertex_prototypes[0]),
        Edge(9, _vertex_prototypes[7], _vertex_prototypes[2]),

        Edge(10, _vertex_prototypes[8], _vertex_prototypes[0]),
        Edge(11, _vertex_prototypes[8], _vertex_prototypes[2]),

    )

    _num_windows = 5
    _num_window_types = 2
