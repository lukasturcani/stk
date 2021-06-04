"""
M6L12 Cube
==========

"""

import numpy as np

from ..cage import Cage
from ..vertices import NonLinearVertex, LinearVertex
from ...topology_graph import Edge


class M6L12Cube(Cage):
    """
    Represents a cage topology graph.

    Metal building blocks with four functional groups are
    required for this topology.

    Ligand building blocks with two functional groups are required for
    this topology.

    When using a :class:`dict` for the `building_blocks` parameter,
    as in :ref:`cage-topology-graph-examples`:
    *Multi-Building Block Cage Construction*, a
    :class:`.BuildingBlock`, with the following number of functional
    groups, needs to be assigned to each of the following vertex ids:

        | 4-functional groups: 0 to 5
        | 2-functional groups: 6 to 17

    See :class:`.Cage` for more details and examples.

    """

    _x = np.sqrt(2)
    _vertex_prototypes = (
        NonLinearVertex(0, [_x, 0, 0]),
        NonLinearVertex(1, [0, _x, 0]),
        NonLinearVertex(2, [-_x, 0, 0]),
        NonLinearVertex(3, [0, -_x, 0]),
        NonLinearVertex(4, [0, 0, _x]),
        NonLinearVertex(5, [0, 0, -_x]),

        LinearVertex(6, [1, 1, 0], False),
        LinearVertex(7, [1, -1, 0], False),
        LinearVertex(8, [1, 0, 1], False),
        LinearVertex(9, [1, 0, -1], False),
        LinearVertex(10, [-1, 1, 0], False),
        LinearVertex(11, [-1, -1, 0], False),
        LinearVertex(12, [-1, 0, 1], False),
        LinearVertex(13, [-1, 0, -1], False),
        LinearVertex(14, [0, 1, 1], False),
        LinearVertex(15, [0, 1, -1], False),
        LinearVertex(16, [0, -1, 1], False),
        LinearVertex(17, [0, -1, -1], False),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[6]),
        Edge(1, _vertex_prototypes[0], _vertex_prototypes[7]),
        Edge(2, _vertex_prototypes[0], _vertex_prototypes[8]),
        Edge(3, _vertex_prototypes[0], _vertex_prototypes[9]),

        Edge(4, _vertex_prototypes[1], _vertex_prototypes[6]),
        Edge(5, _vertex_prototypes[1], _vertex_prototypes[10]),
        Edge(6, _vertex_prototypes[1], _vertex_prototypes[14]),
        Edge(7, _vertex_prototypes[1], _vertex_prototypes[15]),

        Edge(8, _vertex_prototypes[2], _vertex_prototypes[10]),
        Edge(9, _vertex_prototypes[2], _vertex_prototypes[11]),
        Edge(10, _vertex_prototypes[2], _vertex_prototypes[12]),
        Edge(11, _vertex_prototypes[2], _vertex_prototypes[13]),

        Edge(12, _vertex_prototypes[3], _vertex_prototypes[7]),
        Edge(13, _vertex_prototypes[3], _vertex_prototypes[11]),
        Edge(14, _vertex_prototypes[3], _vertex_prototypes[16]),
        Edge(15, _vertex_prototypes[3], _vertex_prototypes[17]),

        Edge(16, _vertex_prototypes[4], _vertex_prototypes[8]),
        Edge(17, _vertex_prototypes[4], _vertex_prototypes[12]),
        Edge(18, _vertex_prototypes[4], _vertex_prototypes[14]),
        Edge(19, _vertex_prototypes[4], _vertex_prototypes[16]),

        Edge(20, _vertex_prototypes[5], _vertex_prototypes[9]),
        Edge(21, _vertex_prototypes[5], _vertex_prototypes[13]),
        Edge(22, _vertex_prototypes[5], _vertex_prototypes[15]),
        Edge(23, _vertex_prototypes[5], _vertex_prototypes[17]),
    )

    _num_windows = 8
    _num_window_types = 1
