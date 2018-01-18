"""
Defines cage topologies from di and 4 functionalised building blocks.

"""


import numpy as np

from .base import _CageTopology,  Vertex, Edge


class TwoPlusFour(_CageTopology):
    """
    Sandwich cage topology from 2 and 4 functionalized building blocks.

    """

    positions_A = [Vertex(0, 0, -10), Vertex(0, 0, 10)]
    alpha, beta = positions_A
    beta.edge_plane_normal = lambda alpha=alpha: np.multiply(
                                        alpha.edge_plane_normal(), -1)

    positions_B = [Edge(alpha, beta),
                   Edge(alpha, beta),
                   Edge(alpha, beta),
                   Edge(alpha, beta)]

    a, b, c, d = positions_B

    a.coord = np.array([10, 0, 0])
    b.coord = np.array([-10, 0, 0])
    c.coord = np.array([0, 10, 0])
    d.coord = np.array([0, -10, 0])

    n_windows = 4
    n_window_types = 1


class ThreePlusSix(_CageTopology):
    """
    A cage topology from 2 and 4 functionalized building blocks.

    """

    positions_A = [Vertex(-20, -10*np.sqrt(3), 0),
                   Vertex(20, -10*np.sqrt(3), 0),
                   Vertex(0, 10*np.sqrt(3), 0)]

    a, b, c = positions_A

    positions_B = [Edge(a, b),
                   Edge(a, b),
                   Edge(b, c),
                   Edge(b, c),
                   Edge(a, c),
                   Edge(a, c)]

    e1, e2, e3, e4, e5, e6 = positions_B
    for e in [e1, e3, e5]:
        e.coord = np.add(e.coord, [0, 0, 10])

    for e in [e2, e4, e6]:
        e.coord = np.add(e.coord, [0, 0, -10])

    n_windows = 5
    n_window_types = 2


class FourPlusEight(_CageTopology):
    """
    A cage topology from 2 and 4 functionalized building blocks.

    """

    positions_A = [Vertex(-10, -10, 0),
                   Vertex(-10, 10, 0),
                   Vertex(10, -10, 0),
                   Vertex(10, 10, 0)]

    a, b, c, d = positions_A

    positions_B = [Edge(a, b),
                   Edge(a, b),
                   Edge(b, d),
                   Edge(b, d),
                   Edge(a, c),
                   Edge(a, c),
                   Edge(c, d),
                   Edge(c, d)]

    e1, e2, e3, e4, e5, e6, e7, e8 = positions_B
    for e in [e1, e3, e5, e7]:
        e.coord = np.add(e.coord, [0, 0, 10])

    for e in [e2, e4, e6, e8]:
        e.coord = np.add(e.coord, [0, 0, -10])

    n_windows = 6
    n_window_types = 2


class SixPlusTwelve(_CageTopology):
    """
    A cage topology from 2 and 4 functionalized building blocks.

    """

    positions_A = [Vertex(-50, -50, 0),
                   Vertex(-50, 50, 0),
                   Vertex(50, -50, 0),
                   Vertex(50, 50, 0),
                   Vertex(0, 0, 50),
                   Vertex(0, 0, -50)]

    a, b, c, d, e, f = positions_A

    positions_B = [Edge(a, b),
                   Edge(b, d),
                   Edge(d, c),
                   Edge(a, c),
                   Edge(e, a),
                   Edge(e, b),
                   Edge(e, c),
                   Edge(e, d),
                   Edge(f, a),
                   Edge(f, b),
                   Edge(f, c),
                   Edge(f, d)]

    n_windows = 8
    n_window_types = 1


class TenPlusTwenty(_CageTopology):
    """
    A cage topology from 2 and 4 functionalized building blocks.

    """

    positions_A = [Vertex(-50, 50, -50),
                   Vertex(-50, -50, -50),
                   Vertex(50, 50, -50),
                   Vertex(50, -50, -50),

                   Vertex(-50, 50, 50),
                   Vertex(-50, -50, 50),
                   Vertex(50, 50, 50),
                   Vertex(50, -50, 50),

                   Vertex(0, 0, 75),
                   Vertex(0, 0, -75)]

    a, b, c, d, e, f, g, h, i, j = positions_A

    positions_B = [Edge(a, c),
                   Edge(a, b),
                   Edge(b, d),
                   Edge(c, d),

                   Edge(e, g),
                   Edge(e, f),
                   Edge(f, h),
                   Edge(g, h),


                   Edge(a, e),
                   Edge(b, f),
                   Edge(c, g),
                   Edge(d, h),

                   Edge(i, e),
                   Edge(i, f),
                   Edge(i, g),
                   Edge(i, h),

                   Edge(j, a),
                   Edge(j, b),
                   Edge(j, c),
                   Edge(j, d)]

    n_windows = 12
    n_window_types = 2
