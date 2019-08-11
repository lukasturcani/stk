"""
Defines cages made from building blocks with 2 and 5 functional groups.

"""

from scipy.constants import golden

from .base import Cage, _CageVertex
from ..topology_graph import Edge


class TwelvePlusThirty(Cage):
    """
    Represents a icosahedron cage topology graph.

    See :class:`.Cage` for more details and examples.

    Attributes
    ----------
    vertices : :class:`tuple` of :class:`.Vertex`
        The vertices which make up the topology graph.

    edges : :class:`tuple` of :class:`.Edge`
        The edges which make up the topology graph.

    """

    # Vertices of a regular origin-centred icosahedron
    # Source: http://eusebeia.dyndns.org/4d/icosahedron
    _vertices = (
        _CageVertex(0, 1, golden),
        _CageVertex(0, -1, golden),
        _CageVertex(0, 1, -golden),
        _CageVertex(0, -1, -golden),
        _CageVertex(1, golden, 0),
        _CageVertex(-1, golden, 0),
        _CageVertex(1, -golden, 0),
        _CageVertex(-1, -golden, 0),
        _CageVertex(golden, 0, 1),
        _CageVertex(-golden, 0, 1),
        _CageVertex(golden, 0, -1),
        _CageVertex(-golden, 0, -1)
    )

    vertices = (
        *_vertices,
        _CageVertex.init_at_center(_vertices[0], _vertices[1]),
        _CageVertex.init_at_center(_vertices[0], _vertices[9]),
        _CageVertex.init_at_center(_vertices[0], _vertices[5]),
        _CageVertex.init_at_center(_vertices[0], _vertices[4]),
        _CageVertex.init_at_center(_vertices[0], _vertices[8]),
        _CageVertex.init_at_center(_vertices[8], _vertices[1]),
        _CageVertex.init_at_center(_vertices[1], _vertices[9]),
        _CageVertex.init_at_center(_vertices[9], _vertices[5]),
        _CageVertex.init_at_center(_vertices[5], _vertices[4]),
        _CageVertex.init_at_center(_vertices[4], _vertices[8]),
        _CageVertex.init_at_center(_vertices[5], _vertices[2]),
        _CageVertex.init_at_center(_vertices[5], _vertices[11]),
        _CageVertex.init_at_center(_vertices[9], _vertices[11]),
        _CageVertex.init_at_center(_vertices[9], _vertices[7]),
        _CageVertex.init_at_center(_vertices[1], _vertices[7]),
        _CageVertex.init_at_center(_vertices[1], _vertices[6]),
        _CageVertex.init_at_center(_vertices[8], _vertices[6]),
        _CageVertex.init_at_center(_vertices[8], _vertices[10]),
        _CageVertex.init_at_center(_vertices[4], _vertices[10]),
        _CageVertex.init_at_center(_vertices[4], _vertices[2]),
        _CageVertex.init_at_center(_vertices[2], _vertices[11]),
        _CageVertex.init_at_center(_vertices[11], _vertices[7]),
        _CageVertex.init_at_center(_vertices[7], _vertices[6]),
        _CageVertex.init_at_center(_vertices[6], _vertices[10]),
        _CageVertex.init_at_center(_vertices[10], _vertices[2]),
        _CageVertex.init_at_center(_vertices[2], _vertices[3]),
        _CageVertex.init_at_center(_vertices[11], _vertices[3]),
        _CageVertex.init_at_center(_vertices[7], _vertices[3]),
        _CageVertex.init_at_center(_vertices[6], _vertices[3]),
        _CageVertex.init_at_center(_vertices[10], _vertices[3])
    )

    edges = (
        Edge(vertices[12], vertices[0]),
        Edge(vertices[12], vertices[1]),
        Edge(vertices[13], vertices[0]),
        Edge(vertices[13], vertices[9]),

        Edge(vertices[14], vertices[0]),
        Edge(vertices[14], vertices[5]),
        Edge(vertices[15], vertices[0]),
        Edge(vertices[15], vertices[4]),

        Edge(vertices[16], vertices[0]),
        Edge(vertices[16], vertices[8]),
        Edge(vertices[17], vertices[8]),
        Edge(vertices[17], vertices[1]),

        Edge(vertices[18], vertices[1]),
        Edge(vertices[18], vertices[9]),
        Edge(vertices[19], vertices[9]),
        Edge(vertices[19], vertices[5]),

        Edge(vertices[20], vertices[5]),
        Edge(vertices[20], vertices[4]),
        Edge(vertices[21], vertices[4]),
        Edge(vertices[21], vertices[8]),

        Edge(vertices[22], vertices[5]),
        Edge(vertices[22], vertices[2]),
        Edge(vertices[23], vertices[5]),
        Edge(vertices[23], vertices[11]),

        Edge(vertices[24], vertices[9]),
        Edge(vertices[24], vertices[11]),
        Edge(vertices[25], vertices[9]),
        Edge(vertices[25], vertices[7]),

        Edge(vertices[26], vertices[1]),
        Edge(vertices[26], vertices[7]),
        Edge(vertices[27], vertices[1]),
        Edge(vertices[27], vertices[6]),

        Edge(vertices[28], vertices[8]),
        Edge(vertices[28], vertices[6]),
        Edge(vertices[29], vertices[8]),
        Edge(vertices[29], vertices[10]),

        Edge(vertices[30], vertices[4]),
        Edge(vertices[30], vertices[10]),
        Edge(vertices[31], vertices[4]),
        Edge(vertices[31], vertices[2]),

        Edge(vertices[32], vertices[2]),
        Edge(vertices[32], vertices[11]),
        Edge(vertices[33], vertices[11]),
        Edge(vertices[33], vertices[7]),

        Edge(vertices[34], vertices[7]),
        Edge(vertices[34], vertices[6]),
        Edge(vertices[35], vertices[6]),
        Edge(vertices[35], vertices[10]),

        Edge(vertices[36], vertices[10]),
        Edge(vertices[36], vertices[2]),
        Edge(vertices[37], vertices[2]),
        Edge(vertices[37], vertices[3]),

        Edge(vertices[38], vertices[11]),
        Edge(vertices[38], vertices[3]),
        Edge(vertices[39], vertices[7]),
        Edge(vertices[39], vertices[3]),

        Edge(vertices[40], vertices[6]),
        Edge(vertices[40], vertices[3]),
        Edge(vertices[41], vertices[10]),
        Edge(vertices[41], vertices[3])
    )

    num_windows = 20
    num_window_types = 1
