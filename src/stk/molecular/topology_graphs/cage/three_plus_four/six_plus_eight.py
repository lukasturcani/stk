"""
Six Plus Eight
==============

"""

from ..cage import Cage
from ..vertices import NonLinearCageVertex
from ...topology_graph import Edge


class SixPlusEight(Cage):
    """
    Represents a cage topology graph.

    Building blocks with three and four functional groups are required
    for this topology graph.

    When using a :class:`dict` for the `building_blocks` parameter,
    as in :ref:`cage-topology-graph-examples`:
    *Multi-Building Block Cage Construction*, a
    :class:`.BuildingBlock`, with the following number of functional
    groups, needs to be assigned to each of the following vertex ids:

        | 4-functional groups: 0 to 5
        | 3-functional groups: 6 to 13

    See :class:`.Cage` for more details and examples.

    """

    _x = 1
    _vertex_prototypes = (
        NonLinearCageVertex(0, [-_x, _x, 0]),
        NonLinearCageVertex(1, [-_x, -_x, 0]),
        NonLinearCageVertex(2, [_x, _x, 0]),
        NonLinearCageVertex(3, [_x, -_x, 0]),

        NonLinearCageVertex(4, [0, 0, _x]),
        NonLinearCageVertex(5, [0, 0, -_x]),
    )

    _vertex_prototypes = (
        *_vertex_prototypes,
        NonLinearCageVertex.init_at_center(
            id=6,
            vertices=(
                _vertex_prototypes[0],
                _vertex_prototypes[4],
                _vertex_prototypes[1]
            ),

        ),
        NonLinearCageVertex.init_at_center(
            id=7,
            vertices=(
                _vertex_prototypes[1],
                _vertex_prototypes[4],
                _vertex_prototypes[3]
            ),
        ),
        NonLinearCageVertex.init_at_center(
            id=8,
            vertices=(
                _vertex_prototypes[4],
                _vertex_prototypes[3],
                _vertex_prototypes[2]
            ),
        ),
        NonLinearCageVertex.init_at_center(
            id=9,
            vertices=(
                _vertex_prototypes[4],
                _vertex_prototypes[2],
                _vertex_prototypes[0]
            ),
        ),

        NonLinearCageVertex.init_at_center(
            id=10,
            vertices=(
                _vertex_prototypes[0],
                _vertex_prototypes[5],
                _vertex_prototypes[1]
            ),
        ),
        NonLinearCageVertex.init_at_center(
            id=11,
            vertices=(
                _vertex_prototypes[5],
                _vertex_prototypes[1],
                _vertex_prototypes[3]
            ),
        ),
        NonLinearCageVertex.init_at_center(
            id=12,
            vertices=(
                _vertex_prototypes[3],
                _vertex_prototypes[5],
                _vertex_prototypes[2]
            ),
        ),
        NonLinearCageVertex.init_at_center(
            id=13,
            vertices=(
                _vertex_prototypes[2],
                _vertex_prototypes[5],
                _vertex_prototypes[0]
            ),
        ),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[6], _vertex_prototypes[0]),
        Edge(1, _vertex_prototypes[6], _vertex_prototypes[4]),
        Edge(2, _vertex_prototypes[6], _vertex_prototypes[1]),

        Edge(3, _vertex_prototypes[7], _vertex_prototypes[1]),
        Edge(4, _vertex_prototypes[7], _vertex_prototypes[4]),
        Edge(5, _vertex_prototypes[7], _vertex_prototypes[3]),

        Edge(6, _vertex_prototypes[8], _vertex_prototypes[4]),
        Edge(7, _vertex_prototypes[8], _vertex_prototypes[3]),
        Edge(8, _vertex_prototypes[8], _vertex_prototypes[2]),

        Edge(9, _vertex_prototypes[9], _vertex_prototypes[4]),
        Edge(10, _vertex_prototypes[9], _vertex_prototypes[2]),
        Edge(11, _vertex_prototypes[9], _vertex_prototypes[0]),

        Edge(12, _vertex_prototypes[10], _vertex_prototypes[0]),
        Edge(13, _vertex_prototypes[10], _vertex_prototypes[5]),
        Edge(14, _vertex_prototypes[10], _vertex_prototypes[1]),

        Edge(15, _vertex_prototypes[11], _vertex_prototypes[5]),
        Edge(16, _vertex_prototypes[11], _vertex_prototypes[1]),
        Edge(17, _vertex_prototypes[11], _vertex_prototypes[3]),

        Edge(18, _vertex_prototypes[12], _vertex_prototypes[3]),
        Edge(19, _vertex_prototypes[12], _vertex_prototypes[5]),
        Edge(20, _vertex_prototypes[12], _vertex_prototypes[2]),

        Edge(21, _vertex_prototypes[13], _vertex_prototypes[2]),
        Edge(22, _vertex_prototypes[13], _vertex_prototypes[5]),
        Edge(23, _vertex_prototypes[13], _vertex_prototypes[0]),

    )

    _num_windows = 12
    _num_window_types = 1
