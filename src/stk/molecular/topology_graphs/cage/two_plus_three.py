"""
Defines cage topologies from 2 and 3 functionalized building blocks.

"""

import numpy as np

from .base import Cage, _CageVertex
from ..topology_graph import Edge


class TwoPlusThree(Cage):
    """
    Represents a cage topology graph.

    See :class:`.Cage` for more details and examples.

    Attributes
    ----------
    vertices : :class:`tuple` of :class:`.Vertex`
        The vertices which make up the topology graph.

    edges : :class:`tuple` of :class:`.Edge`
        The edges which make up the topology graph.

    """

    vertices = (
        _CageVertex(0, 0, 1),
        _CageVertex(0, 0, -1),

        _CageVertex(-1, -0.5*np.sqrt(3), 0, False),
        _CageVertex(1, -0.5*np.sqrt(3), 0, False),
        _CageVertex(0, 0.5*np.sqrt(3), 0, False)

    )

    edges = (
        Edge(vertices[0], vertices[2]),
        Edge(vertices[0], vertices[3]),
        Edge(vertices[0], vertices[4]),
        Edge(vertices[1], vertices[2]),
        Edge(vertices[1], vertices[3]),
        Edge(vertices[1], vertices[4])
    )

    num_windows = 3
    num_window_types = 1


class FourPlusSix(Cage):
    """
    Represents a cage topology graph.

    See :class:`.Cage` for more details and examples.

    Attributes
    ----------
    vertices : :class:`tuple` of :class:`.Vertex`
        The vertices which make up the topology graph.

    edges : :class:`tuple` of :class:`.Edge`
        The edges which make up the topology graph.

    """

    # Vertices of a tetrahdron so that origin is at the origin. Source:
    # http://tinyurl.com/lc262h8.
    _v0, _v1, _v2, _v3 = _vertices = (
        _CageVertex(0, 0, np.sqrt(6)/2),
        _CageVertex(-1, -np.sqrt(3)/3, -np.sqrt(6)/6),
        _CageVertex(1, -np.sqrt(3)/3, -np.sqrt(6)/6),
        _CageVertex(0, 2*np.sqrt(3)/3, -np.sqrt(6)/6)
    )

    vertices = (
        *_vertices,
        _CageVertex.init_at_center(_v0, _v1),
        _CageVertex.init_at_center(_v0, _v2),
        _CageVertex.init_at_center(_v0, _v3),
        _CageVertex.init_at_center(_v1, _v2),
        _CageVertex.init_at_center(_v1, _v3),
        _CageVertex.init_at_center(_v2, _v3)

    )

    edges = (
        Edge(vertices[4], vertices[0]),
        Edge(vertices[4], vertices[1]),
        Edge(vertices[5], vertices[0]),
        Edge(vertices[5], vertices[2]),
        Edge(vertices[6], vertices[0]),
        Edge(vertices[6], vertices[3]),
        Edge(vertices[7], vertices[1]),
        Edge(vertices[7], vertices[2]),
        Edge(vertices[8], vertices[1]),
        Edge(vertices[8], vertices[3]),
        Edge(vertices[9], vertices[2]),
        Edge(vertices[9], vertices[3])
    )

    num_windows = 4
    num_window_types = 1


class FourPlusSix2(Cage):
    """
    Represents a cage topology graph.

    See :class:`.Cage` for more details and examples.

    Attributes
    ----------
    vertices : :class:`tuple` of :class:`.Vertex`
        The vertices which make up the topology graph.

    edges : :class:`tuple` of :class:`.Edge`
        The edges which make up the topology graph.

    """

    _vertices = (
        _CageVertex(1, 0, 1),
        _CageVertex(-1, 0, 1),
        _CageVertex(1, 0, -1),
        _CageVertex(-1, 0, -1),

        _CageVertex(0, -1, 1, False),
        _CageVertex(0, 1, 1, False),
        _CageVertex(0, -1, -1, False),
        _CageVertex(0, 1, -1, False)
    )

    vertices = (
        *_vertices,
        _CageVertex.init_at_center(_vertices[0], _vertices[2]),
        _CageVertex.init_at_center(_vertices[1], _vertices[3])
    )

    edges = (
        Edge(vertices[4], vertices[0]),
        Edge(vertices[4], vertices[1]),
        Edge(vertices[5], vertices[0]),
        Edge(vertices[5], vertices[1]),
        Edge(vertices[6], vertices[2]),
        Edge(vertices[6], vertices[3]),
        Edge(vertices[7], vertices[2]),
        Edge(vertices[7], vertices[3]),
        Edge(vertices[8], vertices[0]),
        Edge(vertices[8], vertices[2]),
        Edge(vertices[9], vertices[1]),
        Edge(vertices[9], vertices[3])
    )

    num_windows = 4
    num_window_types = 2


class SixPlusNine(Cage):
    """
    Represents a cage topology graph.

    See :class:`.Cage` for more details and examples.

    Attributes
    ----------
    vertices : :class:`tuple` of :class:`.Vertex`
        The vertices which make up the topology graph.

    edges : :class:`tuple` of :class:`.Edge`
        The edges which make up the topology graph.

    """

    # source: http://eusebeia.dyndns.org/4d/prism3
    _vertices = (
        _CageVertex(-1, -1/np.sqrt(3), -1),
        _CageVertex(-1, -1/np.sqrt(3), 1),
        _CageVertex(1, -1/np.sqrt(3), -1),
        _CageVertex(1, -1/np.sqrt(3), 1),
        _CageVertex(0, 2/np.sqrt(3), -1),
        _CageVertex(0, 2/np.sqrt(3), 1)
    )
    vertices = (
        *_vertices,
        _CageVertex.init_at_center(_vertices[0], _vertices[1]),
        _CageVertex.init_at_center(_vertices[0], _vertices[2]),
        _CageVertex.init_at_center(_vertices[2], _vertices[3]),
        _CageVertex.init_at_center(_vertices[1], _vertices[3]),
        _CageVertex.init_at_center(_vertices[0], _vertices[4]),
        _CageVertex.init_at_center(_vertices[2], _vertices[4]),
        _CageVertex.init_at_center(_vertices[1], _vertices[5]),
        _CageVertex.init_at_center(_vertices[3], _vertices[5]),
        _CageVertex.init_at_center(_vertices[4], _vertices[5])
    )

    edges = (
        Edge(vertices[6], vertices[0]),
        Edge(vertices[6], vertices[1]),
        Edge(vertices[7], vertices[0]),
        Edge(vertices[7], vertices[2]),
        Edge(vertices[8], vertices[2]),
        Edge(vertices[8], vertices[3]),
        Edge(vertices[9], vertices[1]),
        Edge(vertices[9], vertices[3]),
        Edge(vertices[10], vertices[0]),
        Edge(vertices[10], vertices[4]),
        Edge(vertices[11], vertices[2]),
        Edge(vertices[11], vertices[4]),
        Edge(vertices[12], vertices[1]),
        Edge(vertices[12], vertices[5]),
        Edge(vertices[13], vertices[3]),
        Edge(vertices[13], vertices[5]),
        Edge(vertices[14], vertices[4]),
        Edge(vertices[14], vertices[5])
    )

    num_windows = 5
    num_window_types = 1


class EightPlusTwelve(Cage):
    """
    Represents a cage topology graph.

    See :class:`.Cage` for more details and examples.

    Attributes
    ----------
    vertices : :class:`tuple` of :class:`.Vertex`
        The vertices which make up the topology graph.

    edges : :class:`tuple` of :class:`.Edge`
        The edges which make up the topology graph.

    """

    _vertices = (
        _CageVertex(-1, 1, -1),
        _CageVertex(-1, -1, -1),
        _CageVertex(1, 1, -1),
        _CageVertex(1, -1, -1),

        _CageVertex(-1, 1, 1),
        _CageVertex(-1, -1, 1),
        _CageVertex(1, 1, 1),
        _CageVertex(1, -1, 1)
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
        _CageVertex.init_at_center(_vertices[3], _vertices[7])
    )

    edges = (
        Edge(vertices[8], vertices[0]),
        Edge(vertices[8], vertices[2]),

        Edge(vertices[9], vertices[0]),
        Edge(vertices[9], vertices[1]),

        Edge(vertices[10], vertices[1]),
        Edge(vertices[10], vertices[3]),

        Edge(vertices[11], vertices[2]),
        Edge(vertices[11], vertices[3]),

        Edge(vertices[12], vertices[4]),
        Edge(vertices[12], vertices[6]),

        Edge(vertices[13], vertices[4]),
        Edge(vertices[13], vertices[5]),

        Edge(vertices[14], vertices[5]),
        Edge(vertices[14], vertices[7]),

        Edge(vertices[15], vertices[6]),
        Edge(vertices[15], vertices[7]),

        Edge(vertices[16], vertices[0]),
        Edge(vertices[16], vertices[4]),

        Edge(vertices[17], vertices[1]),
        Edge(vertices[17], vertices[5]),

        Edge(vertices[18], vertices[2]),
        Edge(vertices[18], vertices[6]),

        Edge(vertices[19], vertices[3]),
        Edge(vertices[19], vertices[7])

    )

    num_windows = 6
    num_window_types = 1


class TwentyPlusThirty(Cage):
    """
    Represents a cage topology graph.

    See :class:`.Cage` for more details and examples.

    Attributes
    ----------
    vertices : :class:`tuple` of :class:`.Vertex`
        The vertices which make up the topology graph.

    edges : :class:`tuple` of :class:`.Edge`
        The edges which make up the topology graph.

    """

    # Source: http://tinyurl.com/h2dl949
    _phi = (1 + np.sqrt(5))/2
    _x = 1.5
    _vertices = (
        _CageVertex(_x*_phi, 0.0, _x/_phi),
        _CageVertex(_x*-_phi, 0.0, _x/_phi),
        _CageVertex(_x*-_phi, 0.0, _x/-_phi),
        _CageVertex(_x*_phi, 0.0, _x/-_phi),

        _CageVertex(_x/_phi, _x*_phi, 0.0),
        _CageVertex(_x/_phi, _x*-_phi, 0.0),
        _CageVertex(_x/-_phi, _x*-_phi, 0.0),
        _CageVertex(_x/-_phi, _x*_phi, 0.0),
        _CageVertex(0.0, _x/_phi, _x*_phi),
        _CageVertex(0.0, _x/_phi, _x*-_phi),
        _CageVertex(0.0, _x/-_phi, _x*-_phi),
        _CageVertex(0.0, _x/-_phi, _x*_phi),

        _CageVertex(_x, _x, _x),
        _CageVertex(_x, -_x, _x),
        _CageVertex(-_x, -_x, _x),
        _CageVertex(-_x, _x, _x),
        _CageVertex(-_x, _x, -_x),
        _CageVertex(_x, _x, -_x),
        _CageVertex(_x, -_x, -_x),
        _CageVertex(-_x, -_x, -_x)
    )

    vertices = (
        *_vertices,
        _CageVertex.init_at_center(_vertices[0], _vertices[13]),
        _CageVertex.init_at_center(_vertices[0], _vertices[12]),
        _CageVertex.init_at_center(_vertices[0], _vertices[3]),

        _CageVertex.init_at_center(_vertices[1], _vertices[14]),
        _CageVertex.init_at_center(_vertices[1], _vertices[15]),
        _CageVertex.init_at_center(_vertices[1], _vertices[2]),

        _CageVertex.init_at_center(_vertices[2], _vertices[19]),
        _CageVertex.init_at_center(_vertices[2], _vertices[16]),

        _CageVertex.init_at_center(_vertices[3], _vertices[18]),
        _CageVertex.init_at_center(_vertices[3], _vertices[17]),

        _CageVertex.init_at_center(_vertices[4], _vertices[12]),
        _CageVertex.init_at_center(_vertices[4], _vertices[7]),
        _CageVertex.init_at_center(_vertices[4], _vertices[17]),

        _CageVertex.init_at_center(_vertices[5], _vertices[6]),
        _CageVertex.init_at_center(_vertices[5], _vertices[18]),
        _CageVertex.init_at_center(_vertices[5], _vertices[13]),

        _CageVertex.init_at_center(_vertices[6], _vertices[14]),
        _CageVertex.init_at_center(_vertices[6], _vertices[19]),

        _CageVertex.init_at_center(_vertices[7], _vertices[15]),
        _CageVertex.init_at_center(_vertices[7], _vertices[16]),

        _CageVertex.init_at_center(_vertices[8], _vertices[11]),
        _CageVertex.init_at_center(_vertices[8], _vertices[12]),
        _CageVertex.init_at_center(_vertices[8], _vertices[15]),

        _CageVertex.init_at_center(_vertices[9], _vertices[10]),
        _CageVertex.init_at_center(_vertices[9], _vertices[17]),
        _CageVertex.init_at_center(_vertices[9], _vertices[16]),

        _CageVertex.init_at_center(_vertices[10], _vertices[18]),
        _CageVertex.init_at_center(_vertices[10], _vertices[19]),

        _CageVertex.init_at_center(_vertices[11], _vertices[14]),
        _CageVertex.init_at_center(_vertices[11], _vertices[13])
    )

    edges = (
        Edge(vertices[20], vertices[0]),
        Edge(vertices[20], vertices[13]),
        Edge(vertices[21], vertices[0]),
        Edge(vertices[21], vertices[12]),
        Edge(vertices[22], vertices[0]),
        Edge(vertices[22], vertices[3]),

        Edge(vertices[23], vertices[1]),
        Edge(vertices[23], vertices[14]),
        Edge(vertices[24], vertices[1]),
        Edge(vertices[24], vertices[15]),
        Edge(vertices[25], vertices[1]),
        Edge(vertices[25], vertices[2]),

        Edge(vertices[26], vertices[2]),
        Edge(vertices[26], vertices[19]),
        Edge(vertices[27], vertices[2]),
        Edge(vertices[27], vertices[16]),

        Edge(vertices[28], vertices[3]),
        Edge(vertices[28], vertices[18]),
        Edge(vertices[29], vertices[3]),
        Edge(vertices[29], vertices[17]),

        Edge(vertices[30], vertices[4]),
        Edge(vertices[30], vertices[12]),
        Edge(vertices[31], vertices[4]),
        Edge(vertices[31], vertices[7]),
        Edge(vertices[32], vertices[4]),
        Edge(vertices[32], vertices[17]),

        Edge(vertices[33], vertices[5]),
        Edge(vertices[33], vertices[6]),
        Edge(vertices[34], vertices[5]),
        Edge(vertices[34], vertices[18]),
        Edge(vertices[35], vertices[5]),
        Edge(vertices[35], vertices[13]),

        Edge(vertices[36], vertices[6]),
        Edge(vertices[36], vertices[14]),
        Edge(vertices[37], vertices[6]),
        Edge(vertices[37], vertices[19]),

        Edge(vertices[38], vertices[7]),
        Edge(vertices[38], vertices[15]),
        Edge(vertices[39], vertices[7]),
        Edge(vertices[39], vertices[16]),

        Edge(vertices[40], vertices[8]),
        Edge(vertices[40], vertices[11]),
        Edge(vertices[41], vertices[8]),
        Edge(vertices[41], vertices[12]),
        Edge(vertices[42], vertices[8]),
        Edge(vertices[42], vertices[15]),

        Edge(vertices[43], vertices[9]),
        Edge(vertices[43], vertices[10]),
        Edge(vertices[44], vertices[9]),
        Edge(vertices[44], vertices[17]),
        Edge(vertices[45], vertices[9]),
        Edge(vertices[45], vertices[16]),

        Edge(vertices[46], vertices[10]),
        Edge(vertices[46], vertices[18]),
        Edge(vertices[47], vertices[10]),
        Edge(vertices[47], vertices[19]),

        Edge(vertices[48], vertices[11]),
        Edge(vertices[48], vertices[14]),
        Edge(vertices[49], vertices[11]),
        Edge(vertices[49], vertices[13])
    )

    num_windows = 12
    num_window_types = 1
