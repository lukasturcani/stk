"""
Defines cage topologies from building blocks of 3 and 4 func groups.

"""

from .base import VertexOnlyCageTopology, Vertex


class SixPlusEight(VertexOnlyCageTopology):
    """
    A cage topology of 3 and 4 functional group building blocks.

    """

    x = 1
    positions_A = [Vertex(-x, x, 0),
                   Vertex(-x, -x, 0),
                   Vertex(x, x, 0),
                   Vertex(x, -x, 0),

                   Vertex(0, 0, x),
                   Vertex(0, 0, -x)]

    a, b, c, d, e, f = positions_A

    positions_B = [Vertex.vertex_init(a, e, b),
                   Vertex.vertex_init(b, e, d),
                   Vertex.vertex_init(e, d, c),
                   Vertex.vertex_init(e, c, a),

                   Vertex.vertex_init(a, f, b),
                   Vertex.vertex_init(f, b, d),
                   Vertex.vertex_init(d, f, c),
                   Vertex.vertex_init(c, f, a)]

    n_windows = 12
    n_window_types = 1
