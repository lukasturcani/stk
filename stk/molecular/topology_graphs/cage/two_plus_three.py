"""
Defines cage topologies from 2 and 3 functionalized building blocks.

"""

import numpy as np

from .base import CageTopology, _CageVertex
from ..topology_graph import Edge


class TwoPlusThree(CageTopology):
    """
    A cage topology from 2 and 3 functionalized building blocks.

    """

    vertices = (
        _CageVertex(0, 0, 1),
        _CageVertex(0, 0, -1),

        _CageVertex(-1, -0.5*np.sqrt(3), 0),
        _CageVertex(1, -0.5*np.sqrt(3), 0),
        _CageVertex(0, 0.5*np.sqrt(3), 0)

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


class FourPlusSix(CageTopology):
    """
    Defines the tetrahedral, 4+6, topology.

    This is a topology of cages where 4 building blocks are placed on
    vertices and 6 linkers are placed on the edges between them.

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


class FourPlusSix2(CageTopology):
    """
    Defines the 4+6 topolgy which is not a tetrahedron.

    """

    _vertices = (
        _CageVertex(1, 0, 1),
        _CageVertex(-1, 0, 1),
        _CageVertex(1, 0, -1),
        _CageVertex(-1, 0, -1),

        _CageVertex(0, -1, 1),
        _CageVertex(0, 1, 1),
        _CageVertex(0, -1, -1),
        _CageVertex(0, 1, -1)
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


class SixPlusNine(CageTopology):
    """
    A cage topology from 2 and 3 functionalized building blocks.

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


class EightPlusTwelve(CageTopology):
    """
    A square topology from 2 and 3 functionalized building blocks.

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


class Dodecahedron(CageTopology):
    """
    A dodecahedron cage from 2 and 3 functionalized building blocks.

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
        Edge(_vertices[20], _vertices[0]),
        Edge(_vertices[20], _vertices[13]),
        Edge(_vertices[21], _vertices[0]),
        Edge(_vertices[21], _vertices[12]),
        Edge(_vertices[22], _vertices[0]),
        Edge(_vertices[22], _vertices[3]),

        Edge(_vertices[23], _vertices[1]),
        Edge(_vertices[23], _vertices[14]),
        Edge(_vertices[24], _vertices[1]),
        Edge(_vertices[24], _vertices[15]),
        Edge(_vertices[25], _vertices[1]),
        Edge(_vertices[25], _vertices[2]),

        Edge(_vertices[26], _vertices[2]),
        Edge(_vertices[26], _vertices[19]),
        Edge(_vertices[27], _vertices[2]),
        Edge(_vertices[27], _vertices[16]),

        Edge(_vertices[28], _vertices[3]),
        Edge(_vertices[28], _vertices[18]),
        Edge(_vertices[29], _vertices[3]),
        Edge(_vertices[29], _vertices[17]),

        Edge(_vertices[30], _vertices[4]),
        Edge(_vertices[30], _vertices[12]),
        Edge(_vertices[31], _vertices[4]),
        Edge(_vertices[31], _vertices[7]),
        Edge(_vertices[32], _vertices[4]),
        Edge(_vertices[32], _vertices[17]),

        Edge(_vertices[33], _vertices[5]),
        Edge(_vertices[33], _vertices[6]),
        Edge(_vertices[34], _vertices[5]),
        Edge(_vertices[34], _vertices[18]),
        Edge(_vertices[35], _vertices[5]),
        Edge(_vertices[35], _vertices[13]),

        Edge(_vertices[36], _vertices[6]),
        Edge(_vertices[36], _vertices[14]),
        Edge(_vertices[37], _vertices[6]),
        Edge(_vertices[37], _vertices[19]),

        Edge(_vertices[38], _vertices[7]),
        Edge(_vertices[38], _vertices[15]),
        Edge(_vertices[39], _vertices[7]),
        Edge(_vertices[39], _vertices[16]),

        Edge(_vertices[40], _vertices[8]),
        Edge(_vertices[40], _vertices[11]),
        Edge(_vertices[41], _vertices[8]),
        Edge(_vertices[41], _vertices[12]),
        Edge(_vertices[42], _vertices[8]),
        Edge(_vertices[42], _vertices[15]),

        Edge(_vertices[43], _vertices[9]),
        Edge(_vertices[43], _vertices[10]),
        Edge(_vertices[44], _vertices[9]),
        Edge(_vertices[44], _vertices[17]),
        Edge(_vertices[45], _vertices[9]),
        Edge(_vertices[45], _vertices[16]),

        Edge(_vertices[46], _vertices[10]),
        Edge(_vertices[46], _vertices[18]),
        Edge(_vertices[47], _vertices[10]),
        Edge(_vertices[47], _vertices[19]),

        Edge(_vertices[48], _vertices[11]),
        Edge(_vertices[48], _vertices[14]),
        Edge(_vertices[49], _vertices[11]),
        Edge(_vertices[49], _vertices[13])
    )

    num_windows = 12
    num_window_types = 1
