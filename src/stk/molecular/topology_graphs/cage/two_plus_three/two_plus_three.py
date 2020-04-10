"""
Two Plus Three
==============

"""

import numpy as np

from ..cage import Cage
from ..vertices import _LinearCageVertex, _NonLinearCageVertex
from ...topology_graph import Edge


class TwoPlusThree(Cage):
    """
    Represents a cage topology graph.

    See :class:`.Cage` for more details and examples.

    """

    _vertex_prototypes = (
        _NonLinearCageVertex(0, [0, 0, 1]),
        _NonLinearCageVertex(1, [0, 0, -1]),

        _LinearCageVertex(2, [-1, -0.5*np.sqrt(3), 0], False),
        _LinearCageVertex(3, [1, -0.5*np.sqrt(3), 0], False),
        _LinearCageVertex(4, [0, 0.5*np.sqrt(3), 0], False),

    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[2]),
        Edge(1, _vertex_prototypes[0], _vertex_prototypes[3]),
        Edge(2, _vertex_prototypes[0], _vertex_prototypes[4]),
        Edge(3, _vertex_prototypes[1], _vertex_prototypes[2]),
        Edge(4, _vertex_prototypes[1], _vertex_prototypes[3]),
        Edge(5, _vertex_prototypes[1], _vertex_prototypes[4]),
    )

    _num_windows = 3
    _num_window_types = 1
