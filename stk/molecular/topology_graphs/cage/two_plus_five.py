"""
Defines cage topologies from 2 and 5 functionalized building blocks.

"""

from scipy.constants import golden

from .base import CageTopology, Vertex, Edge


class Icosahedron(CageTopology):
    """
    Defines the icosahederal, 12+30, topology.

    This is a topology of cages where 12 building blocks are placed on
    the vertices and 30 linkers are placed on the edges between them.

    """

    # Vertices of a regular origin-centred icosahedron
    # Source: http://eusebeia.dyndns.org/4d/icosahedron
    positions_A = v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11 = [
                Vertex(0, 1, golden),
                Vertex(0, -1, golden),
                Vertex(0, 1, -golden),
                Vertex(0, -1, -golden),
                Vertex(1, golden, 0),
                Vertex(-1, golden, 0),
                Vertex(1, -golden, 0),
                Vertex(-1, -golden, 0),
                Vertex(golden, 0, 1),
                Vertex(-golden, 0, 1),
                Vertex(golden, 0, -1),
                Vertex(-golden, 0, -1)
    ]

    positions_B = [Edge(v0, v1), Edge(v0, v9),
                   Edge(v0, v5), Edge(v0, v4),
                   Edge(v0, v8), Edge(v8, v1),
                   Edge(v1, v9), Edge(v9, v5),
                   Edge(v5, v4), Edge(v4, v8),
                   Edge(v5, v2), Edge(v5, v11),
                   Edge(v9, v11), Edge(v9, v7),
                   Edge(v1, v7), Edge(v1, v6),
                   Edge(v8, v6), Edge(v8, v10),
                   Edge(v4, v10), Edge(v4, v2),
                   Edge(v2, v11), Edge(v11, v7),
                   Edge(v7, v6), Edge(v6, v10),
                   Edge(v10, v2), Edge(v2, v3),
                   Edge(v11, v3), Edge(v7, v3),
                   Edge(v6, v3), Edge(v10, v3)]

    n_windows = 20
    n_window_types = 1
