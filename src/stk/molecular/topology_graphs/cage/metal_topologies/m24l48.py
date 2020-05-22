"""
M24L48
======

"""

from ..cage import Cage
from ..vertices import _NonLinearCageVertex, _LinearCageVertex
from ...topology_graph import Edge


class M24L48(Cage):
    """
    Represents a cage topology graph.

    Metal building blocks with four functional groups are
    required for this topology.

    Ligand building blocks with two functional groups are required for
    this topology.

    When using a :class:`dict` for initialization, a
    :class:`.BuildingBlock` needs to be assigned to each of the
    following numbers:

        | metals: 0 to 23
        | ligands: 24 to 71

    See :class:`.Cage` for more details and examples.

    """

    _coord1 = 0.621
    _coord2 = -0.621
    _coord3 = _coord1 - (_coord1-_coord2)/2
    _coord4 = -1.5
    _coord5 = 1.5
    _coord6 = _coord5 - (_coord5-_coord1)/2

    _vertex_prototypes = (
        _NonLinearCageVertex(0, [_coord1, _coord1, _coord4]),
        _NonLinearCageVertex(1, [_coord1, _coord1, _coord5]),
        _NonLinearCageVertex(2, [_coord1, _coord4, _coord2]),
        _NonLinearCageVertex(3, [_coord1, _coord5, _coord2]),
        _NonLinearCageVertex(4, [_coord1, _coord2, _coord4]),
        _NonLinearCageVertex(5, [_coord1, _coord2, _coord5]),
        _NonLinearCageVertex(6, [_coord1, _coord4, _coord1]),
        _NonLinearCageVertex(7, [_coord1, _coord5, _coord1]),
        _NonLinearCageVertex(8, [_coord2, _coord1, _coord4]),
        _NonLinearCageVertex(9, [_coord2, _coord1, _coord5]),
        _NonLinearCageVertex(10, [_coord2, _coord2, _coord4]),
        _NonLinearCageVertex(11, [_coord2, _coord2, _coord5]),
        _NonLinearCageVertex(12, [_coord2, _coord4, _coord1]),
        _NonLinearCageVertex(13, [_coord2, _coord5, _coord1]),
        _NonLinearCageVertex(14, [_coord2, _coord4, _coord2]),
        _NonLinearCageVertex(15, [_coord2, _coord5, _coord2]),
        _NonLinearCageVertex(16, [_coord4, _coord1, _coord1]),
        _NonLinearCageVertex(17, [_coord4, _coord2, _coord1]),
        _NonLinearCageVertex(18, [_coord4, _coord2, _coord2]),
        _NonLinearCageVertex(19, [_coord4, _coord1, _coord2]),
        _NonLinearCageVertex(20, [_coord5, _coord1, _coord1]),
        _NonLinearCageVertex(21, [_coord5, _coord2, _coord1]),
        _NonLinearCageVertex(22, [_coord5, _coord2, _coord2]),
        _NonLinearCageVertex(23, [_coord5, _coord1, _coord2]),

        _LinearCageVertex(24, [_coord1, _coord3, _coord4], False),
        _LinearCageVertex(25, [_coord1, _coord4, _coord3], False),
        _LinearCageVertex(26, [_coord1, _coord3, _coord5], False),
        _LinearCageVertex(27, [_coord1, _coord5, _coord3], False),
        _LinearCageVertex(28, [_coord2, _coord3, _coord4], False),
        _LinearCageVertex(29, [_coord2, _coord4, _coord3], False),
        _LinearCageVertex(30, [_coord2, _coord3, _coord5], False),
        _LinearCageVertex(31, [_coord2, _coord5, _coord3], False),

        _LinearCageVertex(32, [_coord3, _coord1, _coord4], False),
        _LinearCageVertex(33, [_coord4, _coord1, _coord3], False),
        _LinearCageVertex(34, [_coord3, _coord1, _coord5], False),
        _LinearCageVertex(35, [_coord5, _coord1, _coord3], False),
        _LinearCageVertex(36, [_coord3, _coord2, _coord4], False),
        _LinearCageVertex(37, [_coord4, _coord2, _coord3], False),
        _LinearCageVertex(38, [_coord3, _coord2, _coord5], False),
        _LinearCageVertex(39, [_coord5, _coord2, _coord3], False),

        _LinearCageVertex(40, [_coord3, _coord4, _coord1], False),
        _LinearCageVertex(41, [_coord4, _coord3, _coord1], False),
        _LinearCageVertex(42, [_coord3, _coord5, _coord1], False),
        _LinearCageVertex(43, [_coord5, _coord3, _coord1], False),
        _LinearCageVertex(44, [_coord3, _coord4, _coord2], False),
        _LinearCageVertex(45, [_coord4, _coord3, _coord2], False),
        _LinearCageVertex(46, [_coord3, _coord5, _coord2], False),
        _LinearCageVertex(47, [_coord5, _coord3, _coord2], False),

        _LinearCageVertex(48, [_coord1, _coord6, _coord6], False),
        _LinearCageVertex(49, [_coord1, _coord6, -_coord6], False),
        _LinearCageVertex(50, [_coord1, -_coord6, _coord6], False),
        _LinearCageVertex(51, [_coord1, -_coord6, -_coord6], False),
        _LinearCageVertex(52, [_coord2, _coord6, _coord6], False),
        _LinearCageVertex(53, [_coord2, _coord6, -_coord6], False),
        _LinearCageVertex(54, [_coord2, -_coord6, _coord6], False),
        _LinearCageVertex(55, [_coord2, -_coord6, -_coord6], False),

        _LinearCageVertex(56, [_coord6, _coord1, _coord6], False),
        _LinearCageVertex(57, [_coord6, _coord1, -_coord6], False),
        _LinearCageVertex(58, [-_coord6, _coord1, _coord6], False),
        _LinearCageVertex(59, [-_coord6, _coord1, -_coord6], False),
        _LinearCageVertex(60, [_coord6, _coord2, _coord6], False),
        _LinearCageVertex(61, [_coord6, _coord2, -_coord6], False),
        _LinearCageVertex(62, [-_coord6, _coord2, _coord6], False),
        _LinearCageVertex(63, [-_coord6, _coord2, -_coord6], False),

        _LinearCageVertex(64, [_coord6, _coord6, _coord1], False),
        _LinearCageVertex(65, [_coord6, -_coord6, _coord1], False),
        _LinearCageVertex(66, [-_coord6, _coord6, _coord1], False),
        _LinearCageVertex(67, [-_coord6, -_coord6, _coord1], False),
        _LinearCageVertex(68, [_coord6, _coord6, _coord2], False),
        _LinearCageVertex(69, [_coord6, -_coord6, _coord2], False),
        _LinearCageVertex(70, [-_coord6, _coord6, _coord2], False),
        _LinearCageVertex(71, [-_coord6, -_coord6, _coord2], False),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[57]),
        Edge(1, _vertex_prototypes[0], _vertex_prototypes[32]),
        Edge(2, _vertex_prototypes[0], _vertex_prototypes[49]),
        Edge(3, _vertex_prototypes[0], _vertex_prototypes[24]),

        Edge(4, _vertex_prototypes[1], _vertex_prototypes[26]),
        Edge(5, _vertex_prototypes[1], _vertex_prototypes[56]),
        Edge(6, _vertex_prototypes[1], _vertex_prototypes[48]),
        Edge(7, _vertex_prototypes[1], _vertex_prototypes[34]),

        Edge(8, _vertex_prototypes[2], _vertex_prototypes[44]),
        Edge(9, _vertex_prototypes[2], _vertex_prototypes[25]),
        Edge(10, _vertex_prototypes[2], _vertex_prototypes[69]),
        Edge(11, _vertex_prototypes[2], _vertex_prototypes[51]),

        Edge(12, _vertex_prototypes[3], _vertex_prototypes[68]),
        Edge(13, _vertex_prototypes[3], _vertex_prototypes[49]),
        Edge(14, _vertex_prototypes[3], _vertex_prototypes[27]),
        Edge(15, _vertex_prototypes[3], _vertex_prototypes[46]),

        Edge(16, _vertex_prototypes[4], _vertex_prototypes[51]),
        Edge(17, _vertex_prototypes[4], _vertex_prototypes[36]),
        Edge(18, _vertex_prototypes[4], _vertex_prototypes[24]),
        Edge(19, _vertex_prototypes[4], _vertex_prototypes[61]),

        Edge(20, _vertex_prototypes[5], _vertex_prototypes[50]),
        Edge(21, _vertex_prototypes[5], _vertex_prototypes[60]),
        Edge(22, _vertex_prototypes[5], _vertex_prototypes[26]),
        Edge(23, _vertex_prototypes[5], _vertex_prototypes[38]),

        Edge(24, _vertex_prototypes[6], _vertex_prototypes[40]),
        Edge(25, _vertex_prototypes[6], _vertex_prototypes[25]),
        Edge(26, _vertex_prototypes[6], _vertex_prototypes[65]),
        Edge(27, _vertex_prototypes[6], _vertex_prototypes[50]),

        Edge(28, _vertex_prototypes[7], _vertex_prototypes[64]),
        Edge(29, _vertex_prototypes[7], _vertex_prototypes[48]),
        Edge(30, _vertex_prototypes[7], _vertex_prototypes[27]),
        Edge(31, _vertex_prototypes[7], _vertex_prototypes[42]),

        Edge(32, _vertex_prototypes[8], _vertex_prototypes[53]),
        Edge(33, _vertex_prototypes[8], _vertex_prototypes[32]),
        Edge(34, _vertex_prototypes[8], _vertex_prototypes[28]),
        Edge(35, _vertex_prototypes[8], _vertex_prototypes[59]),

        Edge(36, _vertex_prototypes[9], _vertex_prototypes[34]),
        Edge(37, _vertex_prototypes[9], _vertex_prototypes[52]),
        Edge(38, _vertex_prototypes[9], _vertex_prototypes[30]),
        Edge(39, _vertex_prototypes[9], _vertex_prototypes[58]),

        Edge(40, _vertex_prototypes[10], _vertex_prototypes[63]),
        Edge(41, _vertex_prototypes[10], _vertex_prototypes[28]),
        Edge(42, _vertex_prototypes[10], _vertex_prototypes[36]),
        Edge(43, _vertex_prototypes[10], _vertex_prototypes[55]),

        Edge(44, _vertex_prototypes[11], _vertex_prototypes[38]),
        Edge(45, _vertex_prototypes[11], _vertex_prototypes[54]),
        Edge(46, _vertex_prototypes[11], _vertex_prototypes[62]),
        Edge(47, _vertex_prototypes[11], _vertex_prototypes[30]),

        Edge(48, _vertex_prototypes[12], _vertex_prototypes[67]),
        Edge(49, _vertex_prototypes[12], _vertex_prototypes[54]),
        Edge(50, _vertex_prototypes[12], _vertex_prototypes[29]),
        Edge(51, _vertex_prototypes[12], _vertex_prototypes[40]),

        Edge(52, _vertex_prototypes[13], _vertex_prototypes[42]),
        Edge(53, _vertex_prototypes[13], _vertex_prototypes[31]),
        Edge(54, _vertex_prototypes[13], _vertex_prototypes[66]),
        Edge(55, _vertex_prototypes[13], _vertex_prototypes[52]),

        Edge(56, _vertex_prototypes[14], _vertex_prototypes[71]),
        Edge(57, _vertex_prototypes[14], _vertex_prototypes[55]),
        Edge(58, _vertex_prototypes[14], _vertex_prototypes[44]),
        Edge(59, _vertex_prototypes[14], _vertex_prototypes[29]),

        Edge(60, _vertex_prototypes[15], _vertex_prototypes[46]),
        Edge(61, _vertex_prototypes[15], _vertex_prototypes[31]),
        Edge(62, _vertex_prototypes[15], _vertex_prototypes[70]),
        Edge(63, _vertex_prototypes[15], _vertex_prototypes[53]),

        Edge(64, _vertex_prototypes[16], _vertex_prototypes[66]),
        Edge(65, _vertex_prototypes[16], _vertex_prototypes[58]),
        Edge(66, _vertex_prototypes[16], _vertex_prototypes[41]),
        Edge(67, _vertex_prototypes[16], _vertex_prototypes[33]),

        Edge(68, _vertex_prototypes[17], _vertex_prototypes[41]),
        Edge(69, _vertex_prototypes[17], _vertex_prototypes[37]),
        Edge(70, _vertex_prototypes[17], _vertex_prototypes[67]),
        Edge(71, _vertex_prototypes[17], _vertex_prototypes[62]),

        Edge(72, _vertex_prototypes[18], _vertex_prototypes[45]),
        Edge(73, _vertex_prototypes[18], _vertex_prototypes[37]),
        Edge(74, _vertex_prototypes[18], _vertex_prototypes[71]),
        Edge(75, _vertex_prototypes[18], _vertex_prototypes[63]),

        Edge(76, _vertex_prototypes[19], _vertex_prototypes[70]),
        Edge(77, _vertex_prototypes[19], _vertex_prototypes[59]),
        Edge(78, _vertex_prototypes[19], _vertex_prototypes[45]),
        Edge(79, _vertex_prototypes[19], _vertex_prototypes[33]),

        Edge(80, _vertex_prototypes[20], _vertex_prototypes[43]),
        Edge(81, _vertex_prototypes[20], _vertex_prototypes[35]),
        Edge(82, _vertex_prototypes[20], _vertex_prototypes[56]),
        Edge(83, _vertex_prototypes[20], _vertex_prototypes[64]),

        Edge(84, _vertex_prototypes[21], _vertex_prototypes[43]),
        Edge(85, _vertex_prototypes[21], _vertex_prototypes[39]),
        Edge(86, _vertex_prototypes[21], _vertex_prototypes[65]),
        Edge(87, _vertex_prototypes[21], _vertex_prototypes[60]),

        Edge(88, _vertex_prototypes[22], _vertex_prototypes[69]),
        Edge(89, _vertex_prototypes[22], _vertex_prototypes[61]),
        Edge(90, _vertex_prototypes[22], _vertex_prototypes[47]),
        Edge(91, _vertex_prototypes[22], _vertex_prototypes[39]),

        Edge(92, _vertex_prototypes[23], _vertex_prototypes[47]),
        Edge(93, _vertex_prototypes[23], _vertex_prototypes[57]),
        Edge(94, _vertex_prototypes[23], _vertex_prototypes[68]),
        Edge(95, _vertex_prototypes[23], _vertex_prototypes[35]),
    )

    _num_windows = 26
    _num_window_types = 2
