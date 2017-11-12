"""
Defines cage topologies from 2 and 3 functionalized building blocks.

"""

import numpy as np

from .base import _CageTopology, Vertex, Edge


class TwoPlusThree(_CageTopology):
    """
    A cage topology from 2 and 3 functionalized building blocks.

    """
    positions_A = [Vertex(0, 0, 20), Vertex(0, 0, -20)]
    a, b = positions_A

    positions_B = [Edge(a, b), Edge(a, b), Edge(a, b)]

    alpha, beta, gamma = positions_B

    b.edge_plane_normal = lambda a=a: np.multiply(
                                            a.edge_plane_normal(), -1)

    alpha.coord = np.array([-20, -10*np.sqrt(3), 0])
    beta.coord = np.array([20, -10*np.sqrt(3), 0])
    gamma.coord = np.array([0, 10*np.sqrt(3), 0])

    n_windows = 3
    n_window_types = 1


class FourPlusSix(_CageTopology):
    """
    Defines the tetrahedral, 4+6, topology.

    This is a topology of cages where 4 building blocks are placed on
    vertices and 6 linkers are placed on the edges between them.

    """

    # Vertices of a tetrahdron so that origin is at the origin. Source:
    # http://tinyurl.com/lc262h8.
    x = 60
    positions_A = v0, v1, v2, v3 = [
                Vertex(0, 0, x*np.sqrt(6)/2),
                Vertex(-x, -x*np.sqrt(3)/3, -x*np.sqrt(6)/6),
                Vertex(x, -x*np.sqrt(3)/3, -x*np.sqrt(6)/6),
                Vertex(0, 2*x*np.sqrt(3)/3, -x*np.sqrt(6)/6)]

    positions_B = [Edge(v0, v1, (0, 1)), Edge(v0, v2, (0, 2)),
                   Edge(v0, v3, (0, 3)), Edge(v1, v2, (1, 2)),
                   Edge(v1, v3, (1, 3)), Edge(v2, v3, (2, 3))]

    n_windows = 4
    n_window_types = 1


class FourPlusSix2(_CageTopology):
    """
    Defines the 4+6 topolgy which is not a tetrahedron.

    """

    positions_A = a, b, c, d = [Vertex(100, 0, 100),
                                Vertex(-100, 0, 100),
                                Vertex(100, 0, -100),
                                Vertex(-100, 0, -100)]

    positions_B = e1, e2, e3, e4, e5, e6 = [Edge(a, b, 'a'),
                                            Edge(a, b, 'b'),
                                            Edge(c, d, 'e'),
                                            Edge(c, d, 'f'),
                                            Edge(a, c, 'd'),
                                            Edge(b, d, 'c')]

    e1.coord = np.array([0, -100, 100])
    e2.coord = np.array([0, 100, 100])
    e3.coord = np.array([0, -100, -100])
    e4.coord = np.array([0, 100, -100])


class SixPlusNine(_CageTopology):
    """
    A cage topology from 2 and 3 functionalized building blocks.

    """

    # source: http://eusebeia.dyndns.org/4d/prism3
    positions_A = [Vertex(-50, -50/np.sqrt(3), -50),
                   Vertex(-50, -50/np.sqrt(3), 50),
                   Vertex(50, -50/np.sqrt(3), -50),
                   Vertex(50, -50/np.sqrt(3), 50),
                   Vertex(0, 100/np.sqrt(3), -50),
                   Vertex(0, 100/np.sqrt(3), 50)]

    a, b, c, d, e, f = positions_A

    positions_B = [Edge(a, b),
                   Edge(a, c),
                   Edge(c, d),
                   Edge(b, d),
                   Edge(a, e),
                   Edge(c, e),
                   Edge(b, f),
                   Edge(d, f),
                   Edge(e, f)]

    n_windows = 5
    n_window_types = 1


class EightPlusTwelve(_CageTopology):
    """
    A square topology from 2 and 3 functionalized building blocks.

    """

    positions_A = [Vertex(-50, 50, -50),
                   Vertex(-50, -50, -50),
                   Vertex(50, 50, -50),
                   Vertex(50, -50, -50),

                   Vertex(-50, 50, 50),
                   Vertex(-50, -50, 50),
                   Vertex(50, 50, 50),
                   Vertex(50, -50, 50)]

    a, b, c, d, e, f, g, h = positions_A

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
                   Edge(d, h)]

    n_windows = 6
    n_window_types = 1


class Dodecahedron(_CageTopology):
    """
    A dodecahedron cage from 2 and 3 functionalized building blocks.

    """

    # Source: http://tinyurl.com/h2dl949
    phi = (1 + np.sqrt(5))/2
    x = 50
    positions_A = [Vertex(x*phi, 0.0, x/phi),
                   Vertex(x*-phi, 0.0, x/phi),
                   Vertex(x*-phi, 0.0, x/-phi),
                   Vertex(x*phi, 0.0, x/-phi),

                   Vertex(x/phi, x*phi, 0.0),
                   Vertex(x/phi, x*-phi, 0.0),
                   Vertex(x/-phi, x*-phi, 0.0),
                   Vertex(x/-phi, x*phi, 0.0),
                   Vertex(0.0, x/phi, x*phi),
                   Vertex(0.0, x/phi, x*-phi),
                   Vertex(0.0, x/-phi, x*-phi),
                   Vertex(0.0, x/-phi, x*phi),

                   Vertex(x, x, x),
                   Vertex(x, -x, x),
                   Vertex(-x, -x, x),
                   Vertex(-x, x, x),
                   Vertex(-x, x, -x),
                   Vertex(x, x, -x),
                   Vertex(x, -x, -x),
                   Vertex(-x, -x, -x)]

    A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T = positions_A
    positions_B = [
         Edge(A, N), Edge(A, M), Edge(A, D), Edge(B, O), Edge(B, P),
         Edge(B, C), Edge(C, T), Edge(C, Q), Edge(D, S), Edge(D, R),
         Edge(E, M), Edge(E, H), Edge(E, R), Edge(F, G), Edge(F, S),
         Edge(F, N), Edge(G, O), Edge(G, T), Edge(H, P), Edge(H, Q),
         Edge(I, L), Edge(I, M), Edge(I, P), Edge(J, K), Edge(J, R),
         Edge(J, Q), Edge(K, S), Edge(K, T), Edge(L, O), Edge(L, N)]

    n_windows = 12
    n_window_types = 1
