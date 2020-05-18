"""
M4L8
====

"""

from ..cage import Cage
from ..vertices import _NonLinearCageVertex, _LinearCageVertex
from ...topology_graph import Edge


class M4L8(Cage):
    """
    Represents a cage topology graph.

    Metal building blocks with four functional groups are
    required for this topology.

    Ligand building blocks with two functional groups are required for
    this topology.

    When using a :class:`dict` for initialization, a
    :class:`.BuildingBlock` needs to be assigned to each of the
    following numbers:

        | metals: (0, 1, 2, 3)
        | ligands: (4, 5, 6, 7, 8, 9, 10, 11)

    See :class:`.Cage` for more details and examples.

    """

    _vertex_prototypes = (
        _NonLinearCageVertex(0, [1, 0, 0]),
        _NonLinearCageVertex(1, [0, 1, 0]),
        _NonLinearCageVertex(2, [-1, 0, 0]),
        _NonLinearCageVertex(3, [0, -1, 0]),

        _LinearCageVertex(4, [1, 1, 1], False),
        _LinearCageVertex(5, [1, 1, -1], False),

        _LinearCageVertex(6, [1, -1, 1], False),
        _LinearCageVertex(7, [1, -1, -1], False),

        _LinearCageVertex(8, [-1, -1, 1], False),
        _LinearCageVertex(9, [-1, -1, -1], False),

        _LinearCageVertex(10, [-1, 1, 1], False),
        _LinearCageVertex(11, [-1, 1, -1], False),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[4]),
        Edge(1, _vertex_prototypes[0], _vertex_prototypes[5]),
        Edge(2, _vertex_prototypes[0], _vertex_prototypes[6]),
        Edge(3, _vertex_prototypes[0], _vertex_prototypes[7]),

        Edge(4, _vertex_prototypes[1], _vertex_prototypes[4]),
        Edge(5, _vertex_prototypes[1], _vertex_prototypes[5]),
        Edge(6, _vertex_prototypes[1], _vertex_prototypes[10]),
        Edge(7, _vertex_prototypes[1], _vertex_prototypes[11]),

        Edge(8, _vertex_prototypes[2], _vertex_prototypes[10]),
        Edge(9, _vertex_prototypes[2], _vertex_prototypes[11]),
        Edge(10, _vertex_prototypes[2], _vertex_prototypes[8]),
        Edge(11, _vertex_prototypes[2], _vertex_prototypes[9]),

        Edge(12, _vertex_prototypes[3], _vertex_prototypes[8]),
        Edge(13, _vertex_prototypes[3], _vertex_prototypes[9]),
        Edge(14, _vertex_prototypes[3], _vertex_prototypes[6]),
        Edge(15, _vertex_prototypes[3], _vertex_prototypes[7]),
    )

    _num_windows = 2
    _num_window_types = 1
