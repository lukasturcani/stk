"""
M8L6 Cube
=========

"""

from ..cage import Cage
from ..vertices import _NonLinearCageVertex
from ...topology_graph import Edge


class M8L6Cube(Cage):
    """
    Represents a cage topology graph.

    Metal building blocks with three functional groups are
    required for this topology.

    Ligand building blocks with four functional groups are required for
    this topology.

    When using a :class:`dict` for initialization, a
    :class:`.BuildingBlock` needs to be assigned to each of the
    following numbers:

        | metals (3 functional groups): 0 to 7
        | ligands (4 functional groups): 8 to 13

    See :class:`.Cage` for more details and examples.

    """

    _vertex_prototypes = (
        _NonLinearCageVertex(0, [2, 2, 2]),
        _NonLinearCageVertex(1, [2, -2, 2]),
        _NonLinearCageVertex(2, [-2, -2, 2]),
        _NonLinearCageVertex(3, [-2, 2, 2]),
        _NonLinearCageVertex(4, [2, 2, -2]),
        _NonLinearCageVertex(5, [2, -2, -2]),
        _NonLinearCageVertex(6, [-2, -2, -2]),
        _NonLinearCageVertex(7, [-2, 2, -2]),
    )

    _vertex_prototypes = (
        *_vertex_prototypes,

        _NonLinearCageVertex(8, position=[0, 0, 2]),
        _NonLinearCageVertex(9, position=[2, 0, 0]),
        _NonLinearCageVertex(10, position=[0, 2, 0]),
        _NonLinearCageVertex(11, position=[-2, 0, 0]),
        _NonLinearCageVertex(12, position=[0, 0, -2]),
        _NonLinearCageVertex(13, position=[0, -2, 0]),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[8]),
        Edge(1, _vertex_prototypes[1], _vertex_prototypes[8]),
        Edge(2, _vertex_prototypes[2], _vertex_prototypes[8]),
        Edge(3, _vertex_prototypes[3], _vertex_prototypes[8]),

        Edge(4, _vertex_prototypes[4], _vertex_prototypes[9]),
        Edge(5, _vertex_prototypes[5], _vertex_prototypes[9]),
        Edge(6, _vertex_prototypes[1], _vertex_prototypes[9]),
        Edge(7, _vertex_prototypes[0], _vertex_prototypes[9]),

        Edge(8, _vertex_prototypes[4], _vertex_prototypes[10]),
        Edge(9, _vertex_prototypes[0], _vertex_prototypes[10]),
        Edge(10, _vertex_prototypes[3], _vertex_prototypes[10]),
        Edge(11, _vertex_prototypes[7], _vertex_prototypes[10]),

        Edge(12, _vertex_prototypes[3], _vertex_prototypes[11]),
        Edge(13, _vertex_prototypes[2], _vertex_prototypes[11]),
        Edge(14, _vertex_prototypes[6], _vertex_prototypes[11]),
        Edge(15, _vertex_prototypes[7], _vertex_prototypes[11]),

        Edge(16, _vertex_prototypes[5], _vertex_prototypes[12]),
        Edge(17, _vertex_prototypes[4], _vertex_prototypes[12]),
        Edge(18, _vertex_prototypes[7], _vertex_prototypes[12]),
        Edge(19, _vertex_prototypes[6], _vertex_prototypes[12]),

        Edge(20, _vertex_prototypes[1], _vertex_prototypes[13]),
        Edge(21, _vertex_prototypes[5], _vertex_prototypes[13]),
        Edge(22, _vertex_prototypes[6], _vertex_prototypes[13]),
        Edge(23, _vertex_prototypes[2], _vertex_prototypes[13]),
    )

    _num_windows = 4
    _num_window_types = 1
