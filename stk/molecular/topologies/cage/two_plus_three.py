"""
Defines cage topologies from 2 and 3 functionalized building blocks.

"""

import numpy as np

from .base import CageTopology, Vertex, Edge


class TwoPlusThree(CageTopology):
    """
    A cage topology from 2 and 3 functionalized building blocks.

    """
    positions_A = a, b = [Vertex(0, 0, 1), Vertex(0, 0, -1)]

    positions_B = alpha, beta, gamma = [Edge(a, b),
                                        Edge(a, b),
                                        Edge(a, b)]

    b.edge_plane_normal = (
        lambda scale, a=a: np.multiply(a.edge_plane_normal(scale), -1)
    )

    alpha.coord = np.array([-1, -0.5*np.sqrt(3), 0])
    beta.coord = np.array([1, -0.5*np.sqrt(3), 0])
    gamma.coord = np.array([0, 0.5*np.sqrt(3), 0])

    n_windows = 3
    n_window_types = 1


class FourPlusSix(CageTopology):
    """
    Defines the tetrahedral, 4+6, topology.

    This is a topology of cages where 4 building blocks are placed on
    vertices and 6 linkers are placed on the edges between them.

    """

    # Vertices of a tetrahdron so that origin is at the origin. Source:
    # http://tinyurl.com/lc262h8.
    positions_A = v0, v1, v2, v3 = [
                Vertex(0, 0, np.sqrt(6)/2),
                Vertex(-1, -np.sqrt(3)/3, -np.sqrt(6)/6),
                Vertex(1, -np.sqrt(3)/3, -np.sqrt(6)/6),
                Vertex(0, 2*np.sqrt(3)/3, -np.sqrt(6)/6)]

    positions_B = [Edge(v0, v1, (0, 1)), Edge(v0, v2, (0, 2)),
                   Edge(v0, v3, (0, 3)), Edge(v1, v2, (1, 2)),
                   Edge(v1, v3, (1, 3)), Edge(v2, v3, (2, 3))]

    n_windows = 4
    n_window_types = 1


class FourPlusSix2(CageTopology):
    """
    Defines the 4+6 topolgy which is not a tetrahedron.

    """

    positions_A = a, b, c, d = [Vertex(1, 0, 1),
                                Vertex(-1, 0, 1),
                                Vertex(1, 0, -1),
                                Vertex(-1, 0, -1)]

    positions_B = e1, e2, e3, e4, e5, e6 = [Edge(a, b, 'a', True),
                                            Edge(a, b, 'b', True),
                                            Edge(c, d, 'e', True),
                                            Edge(c, d, 'f', True),
                                            Edge(a, c, 'd'),
                                            Edge(b, d, 'c')]

    e1.coord = np.array([0, -1, 1])
    e2.coord = np.array([0, 1, 1])
    e3.coord = np.array([0, -1, -1])
    e4.coord = np.array([0, 1, -1])


class SixPlusNine(CageTopology):
    """
    A cage topology from 2 and 3 functionalized building blocks.

    """

    # source: http://eusebeia.dyndns.org/4d/prism3
    positions_A = [Vertex(-1, -1/np.sqrt(3), -1),
                   Vertex(-1, -1/np.sqrt(3), 1),
                   Vertex(1, -1/np.sqrt(3), -1),
                   Vertex(1, -1/np.sqrt(3), 1),
                   Vertex(0, 2/np.sqrt(3), -1),
                   Vertex(0, 2/np.sqrt(3), 1)]

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


class EightPlusTwelve(CageTopology):
    """
    A square topology from 2 and 3 functionalized building blocks.

    """

    positions_A = [Vertex(-1, 1, -1),
                   Vertex(-1, -1, -1),
                   Vertex(1, 1, -1),
                   Vertex(1, -1, -1),

                   Vertex(-1, 1, 1),
                   Vertex(-1, -1, 1),
                   Vertex(1, 1, 1),
                   Vertex(1, -1, 1)]

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


class Dodecahedron(CageTopology):
    """
    A dodecahedron cage from 2 and 3 functionalized building blocks.

    """

    # Source: http://tinyurl.com/h2dl949
    phi = (1 + np.sqrt(5))/2
    x = 1.5
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

    A, B, C, D, E, F, G, *_ = positions_A
    H, I, J, K, L, M, N, O, P, Q, R, S, T = _

    positions_B = [
         Edge(A, N), Edge(A, M), Edge(A, D), Edge(B, O), Edge(B, P),
         Edge(B, C), Edge(C, T), Edge(C, Q), Edge(D, S), Edge(D, R),
         Edge(E, M), Edge(E, H), Edge(E, R), Edge(F, G), Edge(F, S),
         Edge(F, N), Edge(G, O), Edge(G, T), Edge(H, P), Edge(H, Q),
         Edge(I, L), Edge(I, M), Edge(I, P), Edge(J, K), Edge(J, R),
         Edge(J, Q), Edge(K, S), Edge(K, T), Edge(L, O), Edge(L, N)]

    n_windows = 12
    n_window_types = 1
