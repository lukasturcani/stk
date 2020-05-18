"""
M12L24
======

"""

from ..cage import Cage
from ..vertices import _NonLinearCageVertex, _LinearCageVertex
from ...topology_graph import Edge


class M12L24(Cage):
    """
    Represents a cage topology graph.

    Metal building blocks with four functional groups are
    required for this topology.

    Ligand building blocks with two functional groups are required for
    this topology.

    When using a :class:`dict` for initialization, a
    :class:`.BuildingBlock` needs to be assigned to each of the
    following numbers:

        | metals: 0 to 11
        | ligands: 12 to 35

    See :class:`.Cage` for more details and examples.

    """

    _vertex_prototypes = (
        _NonLinearCageVertex(0, [1, 0, 0]),
        _NonLinearCageVertex(1, [-1, 0, 0]),
        _NonLinearCageVertex(2, [0, 1, 0]),
        _NonLinearCageVertex(3, [0, -1, 0]),
        _NonLinearCageVertex(4, [0.5, 0.5, 0.707]),
        _NonLinearCageVertex(5, [0.5, -0.5, 0.707]),
        _NonLinearCageVertex(6, [-0.5, 0.5, 0.707]),
        _NonLinearCageVertex(7, [-0.5, -0.5, 0.707]),
        _NonLinearCageVertex(8, [0.5, 0.5, -0.707]),
        _NonLinearCageVertex(9, [0.5, -0.5, -0.707]),
        _NonLinearCageVertex(10, [-0.5, 0.5, -0.707]),
        _NonLinearCageVertex(11, [-0.5, -0.5, -0.707]),

        _LinearCageVertex(12, [0.9, 0.31, 0.31], False),
        _LinearCageVertex(13, [0.9, 0.31, -0.31], False),
        _LinearCageVertex(14, [0.9, -0.31, 0.31], False),
        _LinearCageVertex(15, [0.9, -0.31, -0.31], False),

        _LinearCageVertex(16, [-0.9, 0.31, 0.31], False),
        _LinearCageVertex(17, [-0.9, 0.31, -0.31], False),
        _LinearCageVertex(18, [-0.9, -0.31, 0.31], False),
        _LinearCageVertex(19, [-0.9, -0.31, -0.31], False),

        _LinearCageVertex(20, [0.31, 0.9, 0.31], False),
        _LinearCageVertex(21, [0.31, 0.9, -0.31], False),
        _LinearCageVertex(22, [-0.31, 0.9, 0.31], False),
        _LinearCageVertex(23, [-0.31, 0.9, -0.31], False),

        _LinearCageVertex(24, [0.31, -0.9, 0.31], False),
        _LinearCageVertex(25, [0.31, -0.9, -0.31], False),
        _LinearCageVertex(26, [-0.31, -0.9, 0.31], False),
        _LinearCageVertex(27, [-0.31, -0.9, -0.31], False),

        _LinearCageVertex(28, [0.58, 0, 0.82], False),
        _LinearCageVertex(29, [-0.58, 0, 0.82], False),
        _LinearCageVertex(30, [0, 0.58, 0.82], False),
        _LinearCageVertex(31, [0, -0.58, 0.82], False),
        _LinearCageVertex(32, [0.58, 0, -0.82], False),
        _LinearCageVertex(33, [-0.58, 0, -0.82], False),
        _LinearCageVertex(34, [0, 0.58, -0.82], False),
        _LinearCageVertex(35, [0, -0.58, -0.82], False),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[12]),
        Edge(1, _vertex_prototypes[0], _vertex_prototypes[13]),
        Edge(2, _vertex_prototypes[0], _vertex_prototypes[14]),
        Edge(3, _vertex_prototypes[0], _vertex_prototypes[15]),

        Edge(4, _vertex_prototypes[1], _vertex_prototypes[16]),
        Edge(5, _vertex_prototypes[1], _vertex_prototypes[17]),
        Edge(6, _vertex_prototypes[1], _vertex_prototypes[18]),
        Edge(7, _vertex_prototypes[1], _vertex_prototypes[19]),

        Edge(8, _vertex_prototypes[2], _vertex_prototypes[20]),
        Edge(9, _vertex_prototypes[2], _vertex_prototypes[21]),
        Edge(10, _vertex_prototypes[2], _vertex_prototypes[22]),
        Edge(11, _vertex_prototypes[2], _vertex_prototypes[23]),

        Edge(12, _vertex_prototypes[3], _vertex_prototypes[24]),
        Edge(13, _vertex_prototypes[3], _vertex_prototypes[25]),
        Edge(14, _vertex_prototypes[3], _vertex_prototypes[26]),
        Edge(15, _vertex_prototypes[3], _vertex_prototypes[27]),

        Edge(16, _vertex_prototypes[4], _vertex_prototypes[28]),
        Edge(17, _vertex_prototypes[4], _vertex_prototypes[30]),
        Edge(18, _vertex_prototypes[4], _vertex_prototypes[12]),
        Edge(19, _vertex_prototypes[4], _vertex_prototypes[20]),

        Edge(20, _vertex_prototypes[5], _vertex_prototypes[14]),
        Edge(21, _vertex_prototypes[5], _vertex_prototypes[24]),
        Edge(22, _vertex_prototypes[5], _vertex_prototypes[28]),
        Edge(23, _vertex_prototypes[5], _vertex_prototypes[31]),

        Edge(24, _vertex_prototypes[6], _vertex_prototypes[16]),
        Edge(25, _vertex_prototypes[6], _vertex_prototypes[29]),
        Edge(26, _vertex_prototypes[6], _vertex_prototypes[30]),
        Edge(27, _vertex_prototypes[6], _vertex_prototypes[22]),

        Edge(28, _vertex_prototypes[7], _vertex_prototypes[18]),
        Edge(29, _vertex_prototypes[7], _vertex_prototypes[26]),
        Edge(30, _vertex_prototypes[7], _vertex_prototypes[31]),
        Edge(31, _vertex_prototypes[7], _vertex_prototypes[29]),

        Edge(32, _vertex_prototypes[8], _vertex_prototypes[13]),
        Edge(33, _vertex_prototypes[8], _vertex_prototypes[32]),
        Edge(34, _vertex_prototypes[8], _vertex_prototypes[34]),
        Edge(35, _vertex_prototypes[8], _vertex_prototypes[21]),

        Edge(36, _vertex_prototypes[9], _vertex_prototypes[15]),
        Edge(37, _vertex_prototypes[9], _vertex_prototypes[32]),
        Edge(38, _vertex_prototypes[9], _vertex_prototypes[35]),
        Edge(39, _vertex_prototypes[9], _vertex_prototypes[25]),

        Edge(40, _vertex_prototypes[10], _vertex_prototypes[17]),
        Edge(41, _vertex_prototypes[10], _vertex_prototypes[23]),
        Edge(42, _vertex_prototypes[10], _vertex_prototypes[34]),
        Edge(43, _vertex_prototypes[10], _vertex_prototypes[33]),

        Edge(44, _vertex_prototypes[11], _vertex_prototypes[19]),
        Edge(45, _vertex_prototypes[11], _vertex_prototypes[33]),
        Edge(46, _vertex_prototypes[11], _vertex_prototypes[27]),
        Edge(47, _vertex_prototypes[11], _vertex_prototypes[35]),
    )

    _num_windows = 14
    _num_window_types = 2
