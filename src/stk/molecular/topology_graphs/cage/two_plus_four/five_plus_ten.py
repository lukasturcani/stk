"""
Five Plus Ten
=============

"""

import numpy as np

from ..cage import Cage
from ..vertices import _LinearCageVertex, _NonLinearCageVertex
from ...topology_graph import Edge


class FivePlusTen(Cage):
    """
    Represents a cage topology graph.

    See :class:`.Cage` for more details and examples.

    """

    _c1 = np.cos(2*np.pi/5)
    _c2 = np.cos(np.pi/5)
    _s1 = np.sin(2*np.pi/5)
    _s2 = np.sin(4*np.pi/5)

    _vertex_prototypes = (
        _NonLinearCageVertex(0, [0, 1, 0], False),
        _NonLinearCageVertex(1, [_s1, _c1, 0], False),
        _NonLinearCageVertex(2, [_s2, -_c2, 0], False),

        _NonLinearCageVertex(3, [-_s2, -_c2, 0], False),
        _NonLinearCageVertex(4, [-_s1, _c1, 0], False),

        _LinearCageVertex(5, [_s1, 1+_c1, 0.5], False),
        _LinearCageVertex(6, [_s1, 1+_c1, -0.5], False),

        _LinearCageVertex(7, [_s1+_s2, _c1-_c2, 0.5], False),
        _LinearCageVertex(8, [_s1+_s2, _c1-_c2, -0.5], False),

        _LinearCageVertex(9, [0, -2*_c2, 0.5], False),
        _LinearCageVertex(10, [0, -2*_c2, -0.5], False),

        _LinearCageVertex(11, [-_s2-_s1, -_c2+_c1, 0.5], False),
        _LinearCageVertex(12, [-_s2-_s1, -_c2+_c1, -0.5], False),

        _LinearCageVertex(13, [-_s1, 1+_c1, 0.5], False),
        _LinearCageVertex(14, [-_s1, 1+_c1, -0.5], False),

    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[5], _vertex_prototypes[0]),
        Edge(1, _vertex_prototypes[5], _vertex_prototypes[1]),
        Edge(2, _vertex_prototypes[6], _vertex_prototypes[0]),
        Edge(3, _vertex_prototypes[6], _vertex_prototypes[1]),

        Edge(4, _vertex_prototypes[7], _vertex_prototypes[1]),
        Edge(5, _vertex_prototypes[7], _vertex_prototypes[2]),
        Edge(6, _vertex_prototypes[8], _vertex_prototypes[1]),
        Edge(7, _vertex_prototypes[8], _vertex_prototypes[2]),

        Edge(8, _vertex_prototypes[9], _vertex_prototypes[2]),
        Edge(9, _vertex_prototypes[9], _vertex_prototypes[3]),
        Edge(10, _vertex_prototypes[10], _vertex_prototypes[2]),
        Edge(11, _vertex_prototypes[10], _vertex_prototypes[3]),

        Edge(12, _vertex_prototypes[11], _vertex_prototypes[3]),
        Edge(13, _vertex_prototypes[11], _vertex_prototypes[4]),
        Edge(14, _vertex_prototypes[12], _vertex_prototypes[3]),
        Edge(15, _vertex_prototypes[12], _vertex_prototypes[4]),

        Edge(16, _vertex_prototypes[13], _vertex_prototypes[4]),
        Edge(17, _vertex_prototypes[13], _vertex_prototypes[0]),
        Edge(18, _vertex_prototypes[14], _vertex_prototypes[4]),
        Edge(19, _vertex_prototypes[14], _vertex_prototypes[0]),

    )

    _num_windows = 7
    _num_window_types = 2
