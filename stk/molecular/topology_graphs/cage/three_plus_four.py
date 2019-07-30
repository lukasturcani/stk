"""
Defines cage topologies from building blocks of 3 and 4 func groups.

"""

from .base import CageTopology, _CageVertex
from ..topology_graph import Edge


class SixPlusEight(CageTopology):
    """
    A cage topology of 3 and 4 functional group building blocks.

    """

    _x = 1
    vertices = (
        _CageVertex(-_x, _x, 0),  # 0
        _CageVertex(-_x, -_x, 0),  # 1
        _CageVertex(_x, _x, 0),  # 2
        _CageVertex(_x, -_x, 0),  # 3
        _CageVertex(0, 0, _x),  # 4
        _CageVertex(0, 0, -_x),  # 5
        _CageVertex(-2*_x, 0, _x),  # 6
        _CageVertex(0, 2*_x, _x),  # 7
        _CageVertex(2*_x, 0, _x),  # 8
        _CageVertex(0, 2*_x, _x),  # 9
        _CageVertex(-2*_x, 0, -_x),  # 10
        _CageVertex(0, -2*_x, -_x),  # 11
        _CageVertex(2*_x, 0, -_x),  # 12
        _CageVertex(0, 2*_x, -_x),  # 13
    )

    edges = (
        Edge(vertices[6], vertices[0]),
        Edge(vertices[6], vertices[5]),
        Edge(vertices[6], vertices[1]),

        Edge(vertices[7], vertices[1]),
        Edge(vertices[7], vertices[5]),
        Edge(vertices[7], vertices[4]),

        Edge(vertices[8], vertices[5]),
        Edge(vertices[8], vertices[4]),
        Edge(vertices[8], vertices[3]),

        Edge(vertices[9], vertices[5]),
        Edge(vertices[9], vertices[3]),
        Edge(vertices[9], vertices[1]),

        Edge(vertices[10], vertices[1]),
        Edge(vertices[10], vertices[2]),
        Edge(vertices[10], vertices[6]),

        Edge(vertices[10], vertices[1]),
        Edge(vertices[10], vertices[2]),
        Edge(vertices[10], vertices[6]),
        Edge(vertices[10], vertices[6]),
    )

    num_windows = 12
    num_window_types = 1
