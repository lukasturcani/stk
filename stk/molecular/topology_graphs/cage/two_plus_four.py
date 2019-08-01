"""
Defines cage topologies from 2 and 4 functionalized building blocks.

"""


import numpy as np

from .base import CageTopology,  _CageVertex
from ..topology_graph import Edge


class TwoPlusFour(CageTopology):
    """
    Sandwich cage topology from 2 and 4 functionalized building blocks.

    """

    vertices = (
        _CageVertex(0, 0, -1),
        _CageVertex(0, 0, 1),

        _CageVertex(2, 0, 0),
        _CageVertex(-2, 0, 0),
        _CageVertex(0, 2, 0),
        _CageVertex(0, -2, 0)
    )

    edges = (
        Edge(vertices[2], vertices[0]),
        Edge(vertices[2], vertices[1]),
        Edge(vertices[3], vertices[0]),
        Edge(vertices[3], vertices[1]),
        Edge(vertices[4], vertices[0]),
        Edge(vertices[4], vertices[1])
    )

    num_windows = 4
    num_window_types = 1


class ThreePlusSix(CageTopology):
    """
    A cage topology from 2 and 4 functionalized building blocks.

    """

    _x = 1
    vertices = (
        _CageVertex(-2*_x, -_x*np.sqrt(3), 0),
        _CageVertex(2*_x, -_x*np.sqrt(3), 0),
        _CageVertex(0, _x*np.sqrt(3), 0),

        _CageVertex(0, -2*_x*np.sqrt(3), _x),
        _CageVertex(0, -2*_x*np.sqrt(3), -_x),

        _CageVertex(2*_x, 0, _x),
        _CageVertex(2*_x, 0, -_x),

        _CageVertex(-2*_x, 0, -_x),
        _CageVertex(-2*_x, 0, -_x),
    )

    edges = (
        Edge(vertices[3], vertices[0]),
        Edge(vertices[3], vertices[1]),

        Edge(vertices[4], vertices[0]),
        Edge(vertices[4], vertices[1]),

        Edge(vertices[5], vertices[1]),
        Edge(vertices[5], vertices[2]),

        Edge(vertices[6], vertices[1]),
        Edge(vertices[6], vertices[2]),

        Edge(vertices[7], vertices[0]),
        Edge(vertices[7], vertices[2]),

        Edge(vertices[8], vertices[0]),
        Edge(vertices[8], vertices[2]),

    )

    num_windows = 5
    num_window_types = 2


class FourPlusEight(CageTopology):
    """
    A cage topology from 2 and 4 functionalized building blocks.

    """

    vertices = (
        _CageVertex(-1, -1, 0),
        _CageVertex(-1, 1, 0),

        _CageVertex(1, -1, 0),
        _CageVertex(1, 1, 0),

        _CageVertex(-2, 0, 1),
        _CageVertex(-2, 0, -1),

        _CageVertex(0, 2, 1),
        _CageVertex(0, 2, -1),

        _CageVertex(0, -2, 1),
        _CageVertex(0, -2, -1),

        _CageVertex(2, 0, 1),
        _CageVertex(2, 0, -1)

    )

    edges = (
        Edge(vertices[4], vertices[0]),
        Edge(vertices[4], vertices[1]),

        Edge(vertices[5], vertices[0]),
        Edge(vertices[5], vertices[1]),

        Edge(vertices[6], vertices[1]),
        Edge(vertices[6], vertices[3]),

        Edge(vertices[7], vertices[1]),
        Edge(vertices[7], vertices[3]),

        Edge(vertices[8], vertices[0]),
        Edge(vertices[8], vertices[2]),

        Edge(vertices[9], vertices[0]),
        Edge(vertices[9], vertices[2]),

        Edge(vertices[10], vertices[2]),
        Edge(vertices[10], vertices[3]),

        Edge(vertices[11], vertices[2]),
        Edge(vertices[11], vertices[3]),

    )

    num_windows = 6
    num_window_types = 2


class FivePlusTen(CageTopology):

    _c1 = np.cos(2*np.pi/5)
    _c2 = np.cos(np.pi/5)
    _s1 = np.sin(2*np.pi/5)
    _s2 = np.sin(4*np.pi/5)

    vertices = (
        _CageVertex(0, 1, 0),
        _CageVertex(_s1, _c1, 0),
        _CageVertex(_s2, -_c2, 0),

        _CageVertex(-_s2, -_c2, 0),
        _CageVertex(-_s1, _c1, 0),

        _CageVertex(_s1, 1+_c1, 0.5),
        _CageVertex(_s1, 1+_c1, -0.5),

        _CageVertex(_s1+_s2, _c1-_c2, 0.5),
        _CageVertex(_s1+_s2, _c1-_c2, -0.5),

        _CageVertex(0, -2*_c2, 0.5),
        _CageVertex(0, -2*_c2, -0.5),

        _CageVertex(-_s2-_s1, -_c2+_c1, 0.5),
        _CageVertex(-_s2-_s1, -_c2+_c1, -0.5),

        _CageVertex(-_s1, 1+_c1, 0.5),
        _CageVertex(-_s1, 1+_c1, -0.5),

    )

    edges = (
        Edge(vertices[5], vertices[0]),
        Edge(vertices[5], vertices[1]),
        Edge(vertices[6], vertices[0]),
        Edge(vertices[6], vertices[1]),

        Edge(vertices[7], vertices[1]),
        Edge(vertices[7], vertices[2]),
        Edge(vertices[8], vertices[1]),
        Edge(vertices[8], vertices[2]),

        Edge(vertices[9], vertices[2]),
        Edge(vertices[9], vertices[3]),
        Edge(vertices[10], vertices[2]),
        Edge(vertices[10], vertices[3]),

        Edge(vertices[11], vertices[3]),
        Edge(vertices[11], vertices[4]),
        Edge(vertices[12], vertices[3]),
        Edge(vertices[12], vertices[4]),

        Edge(vertices[13], vertices[5]),
        Edge(vertices[13], vertices[0]),
        Edge(vertices[14], vertices[5]),
        Edge(vertices[14], vertices[0])

    )

    num_windows = 7
    num_window_types = 2


class SixPlusTwelve(CageTopology):
    """
    A cage topology from 2 and 4 functionalized building blocks.

    """

    _vertices = (
        _CageVertex(-1, -1, 0),
        _CageVertex(-1, 1, 0),
        _CageVertex(1, -1, 0),
        _CageVertex(1, 1, 0),
        _CageVertex(0, 0, 1),
        _CageVertex(0, 0, -1)
    )

    vertices = (
        *_vertices,
        _CageVertex.init_at_center(_vertices[0], _vertices[1]),
        _CageVertex.init_at_center(_vertices[1], _vertices[3]),
        _CageVertex.init_at_center(_vertices[3], _vertices[2]),

        _CageVertex.init_at_center(_vertices[0], _vertices[2]),
        _CageVertex.init_at_center(_vertices[4], _vertices[0]),
        _CageVertex.init_at_center(_vertices[4], _vertices[1]),

        _CageVertex.init_at_center(_vertices[4], _vertices[2]),
        _CageVertex.init_at_center(_vertices[4], _vertices[3]),
        _CageVertex.init_at_center(_vertices[5], _vertices[0]),

        _CageVertex.init_at_center(_vertices[5], _vertices[1]),
        _CageVertex.init_at_center(_vertices[5], _vertices[2]),
        _CageVertex.init_at_center(_vertices[5], _vertices[3])
    )

    edges = (
        Edge(vertices[6], vertices[0]),
        Edge(vertices[6], vertices[1]),
        Edge(vertices[7], vertices[1]),
        Edge(vertices[7], vertices[3]),
        Edge(vertices[8], vertices[3]),
        Edge(vertices[8], vertices[2]),

        Edge(vertices[9], vertices[0]),
        Edge(vertices[9], vertices[2]),
        Edge(vertices[10], vertices[4]),
        Edge(vertices[10], vertices[0]),
        Edge(vertices[11], vertices[4]),
        Edge(vertices[11], vertices[1]),

        Edge(vertices[12], vertices[4]),
        Edge(vertices[12], vertices[2]),
        Edge(vertices[13], vertices[4]),
        Edge(vertices[13], vertices[3]),
        Edge(vertices[14], vertices[5]),
        Edge(vertices[14], vertices[0]),

        Edge(vertices[15], vertices[5]),
        Edge(vertices[15], vertices[1]),
        Edge(vertices[16], vertices[5]),
        Edge(vertices[16], vertices[2]),
        Edge(vertices[17], vertices[5]),
        Edge(vertices[17], vertices[3])
    )

    num_windows = 8
    num_window_types = 1


class EightPlusSixteen(CageTopology):

    _x = 2
    _vertices = (
        _CageVertex(-0.5*_x, 0.5*_x, -0.35*_x),
        _CageVertex(-0.5*_x, -0.5*_x, -0.35*_x),
        _CageVertex(0.5*_x, -0.5*_x, -0.35*_x),
        _CageVertex(0.5*_x, 0.5*_x, -0.35*_x),

        _CageVertex(-_x*np.sqrt(2)/2, 0, _x*0.35),
        _CageVertex(0, -_x*np.sqrt(2)/2, _x*0.35),
        _CageVertex(_x*np.sqrt(2)/2, 0, _x*0.35),
        _CageVertex(0, _x*np.sqrt(2)/2, _x*0.35)
    )

    vertices = (
        *_vertices,

        _CageVertex.init_at_center(_vertices[1], _vertices[5]),
        _CageVertex.init_at_center(_vertices[2], _vertices[5]),
        _CageVertex.init_at_center(_vertices[0], _vertices[4]),
        _CageVertex.init_at_center(_vertices[1], _vertices[4]),

        _CageVertex.init_at_center(_vertices[2], _vertices[6]),
        _CageVertex.init_at_center(_vertices[3], _vertices[6]),
        _CageVertex.init_at_center(_vertices[0], _vertices[7]),
        _CageVertex.init_at_center(_vertices[1], _vertices[7]),

        _CageVertex.init_at_center(_vertices[0], _vertices[1]),
        _CageVertex.init_at_center(_vertices[1], _vertices[2]),
        _CageVertex.init_at_center(_vertices[2], _vertices[3]),
        _CageVertex.init_at_center(_vertices[3], _vertices[0]),

        _CageVertex.init_at_center(_vertices[4], _vertices[5]),
        _CageVertex.init_at_center(_vertices[5], _vertices[6]),
        _CageVertex.init_at_center(_vertices[6], _vertices[7]),
        _CageVertex.init_at_center(_vertices[7], _vertices[4]),

    )

    edges = (
        Edge(vertices[8], vertices[1]),
        Edge(vertices[8], vertices[5]),

        Edge(vertices[9], vertices[2]),
        Edge(vertices[9], vertices[5]),

        Edge(vertices[10], vertices[0]),
        Edge(vertices[10], vertices[4]),

        Edge(vertices[11], vertices[1]),
        Edge(vertices[11], vertices[4]),

        Edge(vertices[12], vertices[2]),
        Edge(vertices[12], vertices[6]),

        Edge(vertices[13], vertices[3]),
        Edge(vertices[13], vertices[6]),

        Edge(vertices[14], vertices[0]),
        Edge(vertices[14], vertices[7]),

        Edge(vertices[15], vertices[1]),
        Edge(vertices[15], vertices[7]),

        Edge(vertices[16], vertices[0]),
        Edge(vertices[16], vertices[1]),

        Edge(vertices[17], vertices[1]),
        Edge(vertices[17], vertices[2]),

        Edge(vertices[18], vertices[2]),
        Edge(vertices[18], vertices[3]),

        Edge(vertices[19], vertices[3]),
        Edge(vertices[19], vertices[0]),

        Edge(vertices[20], vertices[4]),
        Edge(vertices[20], vertices[5]),

        Edge(vertices[21], vertices[5]),
        Edge(vertices[21], vertices[6]),

        Edge(vertices[22], vertices[6]),
        Edge(vertices[22], vertices[7]),

        Edge(vertices[23], vertices[7]),
        Edge(vertices[23], vertices[4])
    )

    num_windows = 10
    num_window_types = 2


class TenPlusTwenty(CageTopology):
    """
    A cage topology from 2 and 4 functionalized building blocks.

    """

    _x = 2
    _vertices = (
        _CageVertex(-_x, _x, -_x),
        _CageVertex(-_x, -_x, -_x),
        _CageVertex(_x, _x, -_x),
        _CageVertex(_x, -_x, -_x),

        _CageVertex(-_x, _x, _x),
        _CageVertex(-_x, -_x, _x),
        _CageVertex(_x, _x, _x),
        _CageVertex(_x, -_x, _x),

        _CageVertex(0, 0, _x*1.5),
        _CageVertex(0, 0, -_x*1.5)
    )

    vertices = (
        *_vertices,

        _CageVertex.init_at_center(_vertices[0], _vertices[2]),
        _CageVertex.init_at_center(_vertices[0], _vertices[1]),
        _CageVertex.init_at_center(_vertices[1], _vertices[3]),
        _CageVertex.init_at_center(_vertices[2], _vertices[3]),

        _CageVertex.init_at_center(_vertices[4], _vertices[6]),
        _CageVertex.init_at_center(_vertices[4], _vertices[5]),
        _CageVertex.init_at_center(_vertices[5], _vertices[7]),
        _CageVertex.init_at_center(_vertices[6], _vertices[7]),

        _CageVertex.init_at_center(_vertices[0], _vertices[4]),
        _CageVertex.init_at_center(_vertices[1], _vertices[5]),
        _CageVertex.init_at_center(_vertices[2], _vertices[6]),
        _CageVertex.init_at_center(_vertices[3], _vertices[7]),

        _CageVertex.init_at_center(_vertices[8], _vertices[4]),
        _CageVertex.init_at_center(_vertices[8], _vertices[5]),
        _CageVertex.init_at_center(_vertices[8], _vertices[6]),
        _CageVertex.init_at_center(_vertices[8], _vertices[7]),

        _CageVertex.init_at_center(_vertices[9], _vertices[0]),
        _CageVertex.init_at_center(_vertices[9], _vertices[1]),
        _CageVertex.init_at_center(_vertices[9], _vertices[2]),
        _CageVertex.init_at_center(_vertices[9], _vertices[3])
    )

    edges = (
        Edge(vertices[10], vertices[0]),
        Edge(vertices[10], vertices[2]),
        Edge(vertices[11], vertices[0]),
        Edge(vertices[11], vertices[1]),

        Edge(vertices[12], vertices[1]),
        Edge(vertices[12], vertices[3]),
        Edge(vertices[13], vertices[2]),
        Edge(vertices[13], vertices[3]),

        Edge(vertices[14], vertices[4]),
        Edge(vertices[14], vertices[6]),
        Edge(vertices[15], vertices[4]),
        Edge(vertices[15], vertices[5]),

        Edge(vertices[16], vertices[5]),
        Edge(vertices[16], vertices[7]),
        Edge(vertices[17], vertices[6]),
        Edge(vertices[17], vertices[7]),

        Edge(vertices[18], vertices[0]),
        Edge(vertices[18], vertices[4]),
        Edge(vertices[19], vertices[1]),
        Edge(vertices[19], vertices[5]),

        Edge(vertices[20], vertices[2]),
        Edge(vertices[20], vertices[6]),
        Edge(vertices[21], vertices[3]),
        Edge(vertices[21], vertices[7]),

        Edge(vertices[22], vertices[8]),
        Edge(vertices[22], vertices[4]),
        Edge(vertices[23], vertices[8]),
        Edge(vertices[23], vertices[5]),

        Edge(vertices[24], vertices[8]),
        Edge(vertices[24], vertices[6]),
        Edge(vertices[25], vertices[8]),
        Edge(vertices[25], vertices[7]),

        Edge(vertices[26], vertices[9]),
        Edge(vertices[26], vertices[0]),
        Edge(vertices[27], vertices[9]),
        Edge(vertices[27], vertices[1]),

        Edge(vertices[28], vertices[9]),
        Edge(vertices[28], vertices[2]),
        Edge(vertices[29], vertices[9]),
        Edge(vertices[29], vertices[3])
    )

    num_windows = 12
    num_window_types = 2
