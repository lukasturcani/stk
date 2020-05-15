"""
M4L4 Tetrahedron
================

"""

import numpy as np

from ..cage import Cage
from ..vertices import _NonLinearCageVertex
from ...topology_graph import Edge


class M4L4Tetrahedron(Cage):
    """
    Represents a cage topology graph.

    Building blocks with three functional groups are required for this
    topology.

    See :class:`.Cage` for more details and examples.

    """

    _vertex_prototypes = (
        _NonLinearCageVertex(0, [0, 0, np.sqrt(6)/2]),
        _NonLinearCageVertex(1, [-1, -np.sqrt(3)/3, -np.sqrt(6)/6]),
        _NonLinearCageVertex(2, [1, -np.sqrt(3)/3, -np.sqrt(6)/6]),
        _NonLinearCageVertex(3, [0, 2*np.sqrt(3)/3, -np.sqrt(6)/6]),
    )

    _vertex_prototypes = (
        *_vertex_prototypes,

        _NonLinearCageVertex.init_at_center(
            id=4,
            vertices=(
                _vertex_prototypes[0],
                _vertex_prototypes[1],
                _vertex_prototypes[2],
            ),
        ),
        _NonLinearCageVertex.init_at_center(
            id=5,
            vertices=(
                _vertex_prototypes[0],
                _vertex_prototypes[1],
                _vertex_prototypes[3],
            ),
        ),
        _NonLinearCageVertex.init_at_center(
            id=6,
            vertices=(
                _vertex_prototypes[0],
                _vertex_prototypes[2],
                _vertex_prototypes[3],
            ),
        ),
        _NonLinearCageVertex.init_at_center(
            id=7,
            vertices=(
                _vertex_prototypes[1],
                _vertex_prototypes[2],
                _vertex_prototypes[3],
            ),
        ),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[4]),
        Edge(1, _vertex_prototypes[0], _vertex_prototypes[5]),
        Edge(2, _vertex_prototypes[0], _vertex_prototypes[6]),
        Edge(3, _vertex_prototypes[1], _vertex_prototypes[4]),
        Edge(4, _vertex_prototypes[1], _vertex_prototypes[5]),
        Edge(5, _vertex_prototypes[1], _vertex_prototypes[7]),
        Edge(6, _vertex_prototypes[2], _vertex_prototypes[4]),
        Edge(7, _vertex_prototypes[2], _vertex_prototypes[6]),
        Edge(8, _vertex_prototypes[2], _vertex_prototypes[7]),
        Edge(9, _vertex_prototypes[3], _vertex_prototypes[5]),
        Edge(10, _vertex_prototypes[3], _vertex_prototypes[6]),
        Edge(11, _vertex_prototypes[3], _vertex_prototypes[7]),
    )

    _num_windows = 4
    _num_window_types = 1
