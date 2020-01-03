"""
Defines cage topologies from 2 and 3 functionalized building blocks.

"""

import numpy as np

from .base import Cage, _CageVertexData
from ..topology_graph import EdgeData


class TwoPlusThree(Cage):
    """
    Represents a cage topology graph.

    See :class:`.Cage` for more details and examples.

    Attributes
    ----------
    vertex_data : :class:`tuple` of :class:`.VertexData`
        A class attribute. Holds the data of the vertices which make up
        the topology graph.

    edge_data : :class:`tuple` of :class:`.EdgeData`
        A class attribute. Holds the data of the edges which make up
        the topology graph.

    vertices : :class:`tuple` of :class:`.Vertex`
        The vertices which make up the topology graph.

    edges : :class:`tuple` of :class:`.Edge`
        The edges which make up the topology graph.

    """

    vertex_data = (
        _CageVertexData(0, 0, 1),
        _CageVertexData(0, 0, -1),

        _CageVertexData(-1, -0.5*np.sqrt(3), 0, False),
        _CageVertexData(1, -0.5*np.sqrt(3), 0, False),
        _CageVertexData(0, 0.5*np.sqrt(3), 0, False)

    )

    edge_data = (
        EdgeData(vertex_data[0], vertex_data[2]),
        EdgeData(vertex_data[0], vertex_data[3]),
        EdgeData(vertex_data[0], vertex_data[4]),
        EdgeData(vertex_data[1], vertex_data[2]),
        EdgeData(vertex_data[1], vertex_data[3]),
        EdgeData(vertex_data[1], vertex_data[4])
    )

    num_windows = 3
    num_window_types = 1


class FourPlusSix(Cage):
    """
    Represents a cage topology graph.

    See :class:`.Cage` for more details and examples.

    Attributes
    ----------
    vertex_data : :class:`tuple` of :class:`.VertexData`
        A class attribute. Holds the data of the vertices which make up
        the topology graph.

    edge_data : :class:`tuple` of :class:`.EdgeData`
        A class attribute. Holds the data of the edges which make up
        the topology graph.

    vertices : :class:`tuple` of :class:`.Vertex`
        The vertices which make up the topology graph.

    edges : :class:`tuple` of :class:`.Edge`
        The edges which make up the topology graph.

    """

    # Vertices of a tetrahdron so that origin is at the origin. Source:
    # http://tinyurl.com/lc262h8.
    _v0, _v1, _v2, _v3 = _vertex_data = (
        _CageVertexData(0, 0, np.sqrt(6)/2),
        _CageVertexData(-1, -np.sqrt(3)/3, -np.sqrt(6)/6),
        _CageVertexData(1, -np.sqrt(3)/3, -np.sqrt(6)/6),
        _CageVertexData(0, 2*np.sqrt(3)/3, -np.sqrt(6)/6)
    )

    vertex_data = (
        *_vertex_data,
        _CageVertexData.init_at_center(_v0, _v1),
        _CageVertexData.init_at_center(_v0, _v2),
        _CageVertexData.init_at_center(_v0, _v3),
        _CageVertexData.init_at_center(_v1, _v2),
        _CageVertexData.init_at_center(_v1, _v3),
        _CageVertexData.init_at_center(_v2, _v3)

    )

    edge_data = (
        EdgeData(vertex_data[4], vertex_data[0]),
        EdgeData(vertex_data[4], vertex_data[1]),
        EdgeData(vertex_data[5], vertex_data[0]),
        EdgeData(vertex_data[5], vertex_data[2]),
        EdgeData(vertex_data[6], vertex_data[0]),
        EdgeData(vertex_data[6], vertex_data[3]),
        EdgeData(vertex_data[7], vertex_data[1]),
        EdgeData(vertex_data[7], vertex_data[2]),
        EdgeData(vertex_data[8], vertex_data[1]),
        EdgeData(vertex_data[8], vertex_data[3]),
        EdgeData(vertex_data[9], vertex_data[2]),
        EdgeData(vertex_data[9], vertex_data[3])
    )

    num_windows = 4
    num_window_types = 1


class FourPlusSix2(Cage):
    """
    Represents a cage topology graph.

    See :class:`.Cage` for more details and examples.

    Attributes
    ----------
    vertex_data : :class:`tuple` of :class:`.VertexData`
        A class attribute. Holds the data of the vertices which make up
        the topology graph.

    edge_data : :class:`tuple` of :class:`.EdgeData`
        A class attribute. Holds the data of the edges which make up
        the topology graph.

    vertices : :class:`tuple` of :class:`.Vertex`
        The vertices which make up the topology graph.

    edges : :class:`tuple` of :class:`.Edge`
        The edges which make up the topology graph.

    """

    _vertex_data = (
        _CageVertexData(1, 0, 1),
        _CageVertexData(-1, 0, 1),
        _CageVertexData(1, 0, -1),
        _CageVertexData(-1, 0, -1),

        _CageVertexData(0, -1, 1, False),
        _CageVertexData(0, 1, 1, False),
        _CageVertexData(0, -1, -1, False),
        _CageVertexData(0, 1, -1, False)
    )

    vertex_data = (
        *_vertex_data,
        _CageVertexData.init_at_center(
            _vertex_data[0], _vertex_data[2]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[1], _vertex_data[3]
        )
    )

    edge_data = (
        EdgeData(vertex_data[4], vertex_data[0]),
        EdgeData(vertex_data[4], vertex_data[1]),
        EdgeData(vertex_data[5], vertex_data[0]),
        EdgeData(vertex_data[5], vertex_data[1]),
        EdgeData(vertex_data[6], vertex_data[2]),
        EdgeData(vertex_data[6], vertex_data[3]),
        EdgeData(vertex_data[7], vertex_data[2]),
        EdgeData(vertex_data[7], vertex_data[3]),
        EdgeData(vertex_data[8], vertex_data[0]),
        EdgeData(vertex_data[8], vertex_data[2]),
        EdgeData(vertex_data[9], vertex_data[1]),
        EdgeData(vertex_data[9], vertex_data[3])
    )

    num_windows = 4
    num_window_types = 2


class SixPlusNine(Cage):
    """
    Represents a cage topology graph.

    See :class:`.Cage` for more details and examples.

    Attributes
    ----------
    vertex_data : :class:`tuple` of :class:`.VertexData`
        A class attribute. Holds the data of the vertices which make up
        the topology graph.

    edge_data : :class:`tuple` of :class:`.EdgeData`
        A class attribute. Holds the data of the edges which make up
        the topology graph.

    vertices : :class:`tuple` of :class:`.Vertex`
        The vertices which make up the topology graph.

    edges : :class:`tuple` of :class:`.Edge`
        The edges which make up the topology graph.

    """

    # source: http://eusebeia.dyndns.org/4d/prism3
    _vertex_data = (
        _CageVertexData(-1, -1/np.sqrt(3), -1),
        _CageVertexData(-1, -1/np.sqrt(3), 1),
        _CageVertexData(1, -1/np.sqrt(3), -1),
        _CageVertexData(1, -1/np.sqrt(3), 1),
        _CageVertexData(0, 2/np.sqrt(3), -1),
        _CageVertexData(0, 2/np.sqrt(3), 1)
    )
    vertex_data = (
        *_vertex_data,
        _CageVertexData.init_at_center(
            _vertex_data[0], _vertex_data[1]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[0], _vertex_data[2]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[2], _vertex_data[3]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[1], _vertex_data[3]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[0], _vertex_data[4]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[2], _vertex_data[4]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[1], _vertex_data[5]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[3], _vertex_data[5]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[4], _vertex_data[5]
        )
    )

    edge_data = (
        EdgeData(vertex_data[6], vertex_data[0]),
        EdgeData(vertex_data[6], vertex_data[1]),
        EdgeData(vertex_data[7], vertex_data[0]),
        EdgeData(vertex_data[7], vertex_data[2]),
        EdgeData(vertex_data[8], vertex_data[2]),
        EdgeData(vertex_data[8], vertex_data[3]),
        EdgeData(vertex_data[9], vertex_data[1]),
        EdgeData(vertex_data[9], vertex_data[3]),
        EdgeData(vertex_data[10], vertex_data[0]),
        EdgeData(vertex_data[10], vertex_data[4]),
        EdgeData(vertex_data[11], vertex_data[2]),
        EdgeData(vertex_data[11], vertex_data[4]),
        EdgeData(vertex_data[12], vertex_data[1]),
        EdgeData(vertex_data[12], vertex_data[5]),
        EdgeData(vertex_data[13], vertex_data[3]),
        EdgeData(vertex_data[13], vertex_data[5]),
        EdgeData(vertex_data[14], vertex_data[4]),
        EdgeData(vertex_data[14], vertex_data[5])
    )

    num_windows = 5
    num_window_types = 1


class EightPlusTwelve(Cage):
    """
    Represents a cage topology graph.

    See :class:`.Cage` for more details and examples.

    Attributes
    ----------
    vertex_data : :class:`tuple` of :class:`.VertexData`
        A class attribute. Holds the data of the vertices which make up
        the topology graph.

    edge_data : :class:`tuple` of :class:`.EdgeData`
        A class attribute. Holds the data of the edges which make up
        the topology graph.

    vertices : :class:`tuple` of :class:`.Vertex`
        The vertices which make up the topology graph.

    edges : :class:`tuple` of :class:`.Edge`
        The edges which make up the topology graph.

    """

    _vertex_data = (
        _CageVertexData(-1, 1, -1),
        _CageVertexData(-1, -1, -1),
        _CageVertexData(1, 1, -1),
        _CageVertexData(1, -1, -1),

        _CageVertexData(-1, 1, 1),
        _CageVertexData(-1, -1, 1),
        _CageVertexData(1, 1, 1),
        _CageVertexData(1, -1, 1)
    )

    vertex_data = (
        *_vertex_data,

        _CageVertexData.init_at_center(
            _vertex_data[0], _vertex_data[2]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[0], _vertex_data[1]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[1], _vertex_data[3]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[2], _vertex_data[3]
        ),

        _CageVertexData.init_at_center(
            _vertex_data[4], _vertex_data[6]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[4], _vertex_data[5]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[5], _vertex_data[7]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[6], _vertex_data[7]
        ),

        _CageVertexData.init_at_center(
            _vertex_data[0], _vertex_data[4]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[1], _vertex_data[5]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[2], _vertex_data[6]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[3], _vertex_data[7]
        )
    )

    edge_data = (
        EdgeData(vertex_data[8], vertex_data[0]),
        EdgeData(vertex_data[8], vertex_data[2]),

        EdgeData(vertex_data[9], vertex_data[0]),
        EdgeData(vertex_data[9], vertex_data[1]),

        EdgeData(vertex_data[10], vertex_data[1]),
        EdgeData(vertex_data[10], vertex_data[3]),

        EdgeData(vertex_data[11], vertex_data[2]),
        EdgeData(vertex_data[11], vertex_data[3]),

        EdgeData(vertex_data[12], vertex_data[4]),
        EdgeData(vertex_data[12], vertex_data[6]),

        EdgeData(vertex_data[13], vertex_data[4]),
        EdgeData(vertex_data[13], vertex_data[5]),

        EdgeData(vertex_data[14], vertex_data[5]),
        EdgeData(vertex_data[14], vertex_data[7]),

        EdgeData(vertex_data[15], vertex_data[6]),
        EdgeData(vertex_data[15], vertex_data[7]),

        EdgeData(vertex_data[16], vertex_data[0]),
        EdgeData(vertex_data[16], vertex_data[4]),

        EdgeData(vertex_data[17], vertex_data[1]),
        EdgeData(vertex_data[17], vertex_data[5]),

        EdgeData(vertex_data[18], vertex_data[2]),
        EdgeData(vertex_data[18], vertex_data[6]),

        EdgeData(vertex_data[19], vertex_data[3]),
        EdgeData(vertex_data[19], vertex_data[7])

    )

    num_windows = 6
    num_window_types = 1


class TwentyPlusThirty(Cage):
    """
    Represents a cage topology graph.

    See :class:`.Cage` for more details and examples.

    Attributes
    ----------
    vertex_data : :class:`tuple` of :class:`.VertexData`
        A class attribute. Holds the data of the vertices which make up
        the topology graph.

    edge_data : :class:`tuple` of :class:`.EdgeData`
        A class attribute. Holds the data of the edges which make up
        the topology graph.

    vertices : :class:`tuple` of :class:`.Vertex`
        The vertices which make up the topology graph.

    edges : :class:`tuple` of :class:`.Edge`
        The edges which make up the topology graph.

    """

    # Source: http://tinyurl.com/h2dl949
    _phi = (1 + np.sqrt(5))/2
    _x = 1.5
    _vertex_data = (
        _CageVertexData(_x*_phi, 0.0, _x/_phi),
        _CageVertexData(_x*-_phi, 0.0, _x/_phi),
        _CageVertexData(_x*-_phi, 0.0, _x/-_phi),
        _CageVertexData(_x*_phi, 0.0, _x/-_phi),

        _CageVertexData(_x/_phi, _x*_phi, 0.0),
        _CageVertexData(_x/_phi, _x*-_phi, 0.0),
        _CageVertexData(_x/-_phi, _x*-_phi, 0.0),
        _CageVertexData(_x/-_phi, _x*_phi, 0.0),
        _CageVertexData(0.0, _x/_phi, _x*_phi),
        _CageVertexData(0.0, _x/_phi, _x*-_phi),
        _CageVertexData(0.0, _x/-_phi, _x*-_phi),
        _CageVertexData(0.0, _x/-_phi, _x*_phi),

        _CageVertexData(_x, _x, _x),
        _CageVertexData(_x, -_x, _x),
        _CageVertexData(-_x, -_x, _x),
        _CageVertexData(-_x, _x, _x),
        _CageVertexData(-_x, _x, -_x),
        _CageVertexData(_x, _x, -_x),
        _CageVertexData(_x, -_x, -_x),
        _CageVertexData(-_x, -_x, -_x)
    )

    vertex_data = (
        *_vertex_data,
        _CageVertexData.init_at_center(
            _vertex_data[0], _vertex_data[13]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[0], _vertex_data[12]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[0], _vertex_data[3]
        ),

        _CageVertexData.init_at_center(
            _vertex_data[1], _vertex_data[14]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[1], _vertex_data[15]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[1], _vertex_data[2]
        ),

        _CageVertexData.init_at_center(
            _vertex_data[2], _vertex_data[19]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[2], _vertex_data[16]
        ),

        _CageVertexData.init_at_center(
            _vertex_data[3], _vertex_data[18]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[3], _vertex_data[17]
        ),

        _CageVertexData.init_at_center(
            _vertex_data[4], _vertex_data[12]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[4], _vertex_data[7]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[4], _vertex_data[17]
        ),

        _CageVertexData.init_at_center(
            _vertex_data[5], _vertex_data[6]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[5], _vertex_data[18]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[5], _vertex_data[13]
        ),

        _CageVertexData.init_at_center(
            _vertex_data[6], _vertex_data[14]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[6], _vertex_data[19]
        ),

        _CageVertexData.init_at_center(
            _vertex_data[7], _vertex_data[15]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[7], _vertex_data[16]
        ),

        _CageVertexData.init_at_center(
            _vertex_data[8], _vertex_data[11]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[8], _vertex_data[12]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[8], _vertex_data[15]
        ),

        _CageVertexData.init_at_center(
            _vertex_data[9], _vertex_data[10]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[9], _vertex_data[17]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[9], _vertex_data[16]
        ),

        _CageVertexData.init_at_center(
            _vertex_data[10], _vertex_data[18]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[10], _vertex_data[19]
        ),

        _CageVertexData.init_at_center(
            _vertex_data[11], _vertex_data[14]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[11], _vertex_data[13])
    )

    edge_data = (
        EdgeData(vertex_data[20], vertex_data[0]),
        EdgeData(vertex_data[20], vertex_data[13]),
        EdgeData(vertex_data[21], vertex_data[0]),
        EdgeData(vertex_data[21], vertex_data[12]),
        EdgeData(vertex_data[22], vertex_data[0]),
        EdgeData(vertex_data[22], vertex_data[3]),

        EdgeData(vertex_data[23], vertex_data[1]),
        EdgeData(vertex_data[23], vertex_data[14]),
        EdgeData(vertex_data[24], vertex_data[1]),
        EdgeData(vertex_data[24], vertex_data[15]),
        EdgeData(vertex_data[25], vertex_data[1]),
        EdgeData(vertex_data[25], vertex_data[2]),

        EdgeData(vertex_data[26], vertex_data[2]),
        EdgeData(vertex_data[26], vertex_data[19]),
        EdgeData(vertex_data[27], vertex_data[2]),
        EdgeData(vertex_data[27], vertex_data[16]),

        EdgeData(vertex_data[28], vertex_data[3]),
        EdgeData(vertex_data[28], vertex_data[18]),
        EdgeData(vertex_data[29], vertex_data[3]),
        EdgeData(vertex_data[29], vertex_data[17]),

        EdgeData(vertex_data[30], vertex_data[4]),
        EdgeData(vertex_data[30], vertex_data[12]),
        EdgeData(vertex_data[31], vertex_data[4]),
        EdgeData(vertex_data[31], vertex_data[7]),
        EdgeData(vertex_data[32], vertex_data[4]),
        EdgeData(vertex_data[32], vertex_data[17]),

        EdgeData(vertex_data[33], vertex_data[5]),
        EdgeData(vertex_data[33], vertex_data[6]),
        EdgeData(vertex_data[34], vertex_data[5]),
        EdgeData(vertex_data[34], vertex_data[18]),
        EdgeData(vertex_data[35], vertex_data[5]),
        EdgeData(vertex_data[35], vertex_data[13]),

        EdgeData(vertex_data[36], vertex_data[6]),
        EdgeData(vertex_data[36], vertex_data[14]),
        EdgeData(vertex_data[37], vertex_data[6]),
        EdgeData(vertex_data[37], vertex_data[19]),

        EdgeData(vertex_data[38], vertex_data[7]),
        EdgeData(vertex_data[38], vertex_data[15]),
        EdgeData(vertex_data[39], vertex_data[7]),
        EdgeData(vertex_data[39], vertex_data[16]),

        EdgeData(vertex_data[40], vertex_data[8]),
        EdgeData(vertex_data[40], vertex_data[11]),
        EdgeData(vertex_data[41], vertex_data[8]),
        EdgeData(vertex_data[41], vertex_data[12]),
        EdgeData(vertex_data[42], vertex_data[8]),
        EdgeData(vertex_data[42], vertex_data[15]),

        EdgeData(vertex_data[43], vertex_data[9]),
        EdgeData(vertex_data[43], vertex_data[10]),
        EdgeData(vertex_data[44], vertex_data[9]),
        EdgeData(vertex_data[44], vertex_data[17]),
        EdgeData(vertex_data[45], vertex_data[9]),
        EdgeData(vertex_data[45], vertex_data[16]),

        EdgeData(vertex_data[46], vertex_data[10]),
        EdgeData(vertex_data[46], vertex_data[18]),
        EdgeData(vertex_data[47], vertex_data[10]),
        EdgeData(vertex_data[47], vertex_data[19]),

        EdgeData(vertex_data[48], vertex_data[11]),
        EdgeData(vertex_data[48], vertex_data[14]),
        EdgeData(vertex_data[49], vertex_data[11]),
        EdgeData(vertex_data[49], vertex_data[13])
    )

    num_windows = 12
    num_window_types = 1
