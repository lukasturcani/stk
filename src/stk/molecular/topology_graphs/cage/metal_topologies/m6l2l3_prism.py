"""
M6L2L3 Prism
============

"""

import numpy as np

from ..cage import Cage
from ..vertices import _NonLinearCageVertex
from ...topology_graph import Edge


class M6L2L3Prism(Cage):
    """
    Represents a cage topology graph.

    Metal building blocks with three functional groups are
    required for this topology.

    Ligand building blocks with three and four functional groups are
    required for ligand type A and B, respectively, on this topology.

    When using a :class:`dict` for initialization, a
    :class:`.BuildingBlock` needs to be assigned to each of the
    following numbers:

        | metals: 0 to 5
        | ligand A (3 functional groups): (6, 7)
        | ligand B (4 functional groups): (8, 9, 10)

    See :class:`.Cage` for more details and examples.

    """

    _vertex_prototypes = (
        _NonLinearCageVertex(0, [-1, -1/np.sqrt(3), 1]),
        _NonLinearCageVertex(1, [1, -1/np.sqrt(3), 1]),
        _NonLinearCageVertex(2, [0, 2/np.sqrt(3), 1]),

        _NonLinearCageVertex(3, [-1, -1/np.sqrt(3), -1]),
        _NonLinearCageVertex(4, [1, -1/np.sqrt(3), -1]),
        _NonLinearCageVertex(5, [0, 2/np.sqrt(3), -1]),
    )

    _vertex_prototypes = (
        *_vertex_prototypes,

        _NonLinearCageVertex.init_at_center(
            id=6,
            vertices=(
                _vertex_prototypes[0],
                _vertex_prototypes[1],
                _vertex_prototypes[2],
            ),
        ),
        _NonLinearCageVertex.init_at_center(
            id=7,
            vertices=(
                _vertex_prototypes[3],
                _vertex_prototypes[4],
                _vertex_prototypes[5],
            ),
        ),

        _NonLinearCageVertex.init_at_center(
            id=8,
            vertices=(
                _vertex_prototypes[0],
                _vertex_prototypes[1],
                _vertex_prototypes[3],
                _vertex_prototypes[4],
            ),
        ),
        _NonLinearCageVertex.init_at_center(
            id=9,
            vertices=(
                _vertex_prototypes[1],
                _vertex_prototypes[2],
                _vertex_prototypes[4],
                _vertex_prototypes[5],
            ),
        ),
        _NonLinearCageVertex.init_at_center(
            id=10,
            vertices=(
                _vertex_prototypes[2],
                _vertex_prototypes[0],
                _vertex_prototypes[5],
                _vertex_prototypes[3],
            ),
        ),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[6]),
        Edge(1, _vertex_prototypes[0], _vertex_prototypes[8]),
        Edge(2, _vertex_prototypes[0], _vertex_prototypes[10]),
        Edge(3, _vertex_prototypes[1], _vertex_prototypes[6]),
        Edge(4, _vertex_prototypes[1], _vertex_prototypes[8]),
        Edge(5, _vertex_prototypes[1], _vertex_prototypes[9]),
        Edge(6, _vertex_prototypes[2], _vertex_prototypes[6]),
        Edge(7, _vertex_prototypes[2], _vertex_prototypes[9]),
        Edge(8, _vertex_prototypes[2], _vertex_prototypes[10]),
        Edge(9, _vertex_prototypes[3], _vertex_prototypes[7]),
        Edge(10, _vertex_prototypes[3], _vertex_prototypes[8]),
        Edge(11, _vertex_prototypes[3], _vertex_prototypes[10]),
        Edge(12, _vertex_prototypes[4], _vertex_prototypes[7]),
        Edge(13, _vertex_prototypes[4], _vertex_prototypes[8]),
        Edge(14, _vertex_prototypes[4], _vertex_prototypes[9]),
        Edge(15, _vertex_prototypes[5], _vertex_prototypes[7]),
        Edge(16, _vertex_prototypes[5], _vertex_prototypes[9]),
        Edge(17, _vertex_prototypes[5], _vertex_prototypes[10]),
    )

    _num_windows = 4
    _num_window_types = 1
