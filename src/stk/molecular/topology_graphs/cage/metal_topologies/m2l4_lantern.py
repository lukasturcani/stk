"""
M2L4 Lantern
============

"""

from ..cage import Cage
from ..vertices import _NonLinearCageVertex, _LinearCageVertex
from ...topology_graph import Edge


class M2L4Lantern(Cage):
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

        | 4-functional groups: 0 to 1
        | 2-functional groups: 2 to 5

    See :class:`.Cage` for more details and examples.

    """

    _vertex_prototypes = (
        _NonLinearCageVertex(0, [0, 0.5, 0]),
        _NonLinearCageVertex(1, [0, -0.5, 0]),

        _LinearCageVertex(2, [1, 0, 0], False),
        _LinearCageVertex(3, [0, 0, 1], False),
        _LinearCageVertex(4, [-1, 0, 0], False),
        _LinearCageVertex(5, [0, 0, -1], False),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[2]),
        Edge(1, _vertex_prototypes[0], _vertex_prototypes[3]),
        Edge(2, _vertex_prototypes[0], _vertex_prototypes[4]),
        Edge(3, _vertex_prototypes[0], _vertex_prototypes[5]),

        Edge(4, _vertex_prototypes[1], _vertex_prototypes[2]),
        Edge(5, _vertex_prototypes[1], _vertex_prototypes[3]),
        Edge(6, _vertex_prototypes[1], _vertex_prototypes[4]),
        Edge(7, _vertex_prototypes[1], _vertex_prototypes[5]),
    )

    _num_windows = 4
    _num_window_types = 1
