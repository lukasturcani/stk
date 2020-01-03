"""
Defines cages made from building blocks with 2 and 4 functional groups.

"""


import numpy as np

from .base import Cage,  _CageVertexData
from ..topology_graph import EdgeData


class TwoPlusFour(Cage):
    """
    Represents a capsule cage topology graph.

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
        _CageVertexData(0, 0, -1),
        _CageVertexData(0, 0, 1),

        _CageVertexData(2, 0, 0, False),
        _CageVertexData(-2, 0, 0, False),
        _CageVertexData(0, 2, 0, False),
        _CageVertexData(0, -2, 0, False)
    )

    edge_data = (
        EdgeData(vertex_data[2], vertex_data[0]),
        EdgeData(vertex_data[2], vertex_data[1]),
        EdgeData(vertex_data[3], vertex_data[0]),
        EdgeData(vertex_data[3], vertex_data[1]),
        EdgeData(vertex_data[4], vertex_data[0]),
        EdgeData(vertex_data[4], vertex_data[1]),
        EdgeData(vertex_data[5], vertex_data[0]),
        EdgeData(vertex_data[5], vertex_data[1])
    )

    num_windows = 4
    num_window_types = 1


class ThreePlusSix(Cage):
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

    _x = 1
    vertex_data = (
        _CageVertexData(-2*_x, -_x*np.sqrt(3), 0, False),
        _CageVertexData(2*_x, -_x*np.sqrt(3), 0, False),
        _CageVertexData(0, _x*np.sqrt(3), 0, False),

        _CageVertexData(0, -2*_x*np.sqrt(3), _x, False),
        _CageVertexData(0, -2*_x*np.sqrt(3), -_x, False),

        _CageVertexData(2*_x, 0, _x, False),
        _CageVertexData(2*_x, 0, -_x, False),

        _CageVertexData(-2*_x, 0, _x, False),
        _CageVertexData(-2*_x, 0, -_x, False),
    )

    edge_data = (
        EdgeData(vertex_data[3], vertex_data[0]),
        EdgeData(vertex_data[3], vertex_data[1]),

        EdgeData(vertex_data[4], vertex_data[0]),
        EdgeData(vertex_data[4], vertex_data[1]),

        EdgeData(vertex_data[5], vertex_data[1]),
        EdgeData(vertex_data[5], vertex_data[2]),

        EdgeData(vertex_data[6], vertex_data[1]),
        EdgeData(vertex_data[6], vertex_data[2]),

        EdgeData(vertex_data[7], vertex_data[0]),
        EdgeData(vertex_data[7], vertex_data[2]),

        EdgeData(vertex_data[8], vertex_data[0]),
        EdgeData(vertex_data[8], vertex_data[2]),

    )

    num_windows = 5
    num_window_types = 2


class FourPlusEight(Cage):
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
        _CageVertexData(-1, -1, 0, False),
        _CageVertexData(-1, 1, 0, False),

        _CageVertexData(1, -1, 0, False),
        _CageVertexData(1, 1, 0, False),

        _CageVertexData(-2, 0, 1, False),
        _CageVertexData(-2, 0, -1, False),

        _CageVertexData(0, 2, 1, False),
        _CageVertexData(0, 2, -1, False),

        _CageVertexData(0, -2, 1, False),
        _CageVertexData(0, -2, -1, False),

        _CageVertexData(2, 0, 1, False),
        _CageVertexData(2, 0, -1, False)

    )

    edge_data = (
        EdgeData(vertex_data[4], vertex_data[0]),
        EdgeData(vertex_data[4], vertex_data[1]),

        EdgeData(vertex_data[5], vertex_data[0]),
        EdgeData(vertex_data[5], vertex_data[1]),

        EdgeData(vertex_data[6], vertex_data[1]),
        EdgeData(vertex_data[6], vertex_data[3]),

        EdgeData(vertex_data[7], vertex_data[1]),
        EdgeData(vertex_data[7], vertex_data[3]),

        EdgeData(vertex_data[8], vertex_data[0]),
        EdgeData(vertex_data[8], vertex_data[2]),

        EdgeData(vertex_data[9], vertex_data[0]),
        EdgeData(vertex_data[9], vertex_data[2]),

        EdgeData(vertex_data[10], vertex_data[2]),
        EdgeData(vertex_data[10], vertex_data[3]),

        EdgeData(vertex_data[11], vertex_data[2]),
        EdgeData(vertex_data[11], vertex_data[3]),

    )

    num_windows = 6
    num_window_types = 2


class FivePlusTen(Cage):
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

    _c1 = np.cos(2*np.pi/5)
    _c2 = np.cos(np.pi/5)
    _s1 = np.sin(2*np.pi/5)
    _s2 = np.sin(4*np.pi/5)

    vertex_data = (
        _CageVertexData(0, 1, 0, False),
        _CageVertexData(_s1, _c1, 0, False),
        _CageVertexData(_s2, -_c2, 0, False),

        _CageVertexData(-_s2, -_c2, 0, False),
        _CageVertexData(-_s1, _c1, 0, False),

        _CageVertexData(_s1, 1+_c1, 0.5, False),
        _CageVertexData(_s1, 1+_c1, -0.5, False),

        _CageVertexData(_s1+_s2, _c1-_c2, 0.5, False),
        _CageVertexData(_s1+_s2, _c1-_c2, -0.5, False),

        _CageVertexData(0, -2*_c2, 0.5, False),
        _CageVertexData(0, -2*_c2, -0.5, False),

        _CageVertexData(-_s2-_s1, -_c2+_c1, 0.5, False),
        _CageVertexData(-_s2-_s1, -_c2+_c1, -0.5, False),

        _CageVertexData(-_s1, 1+_c1, 0.5, False),
        _CageVertexData(-_s1, 1+_c1, -0.5, False),

    )

    edge_data = (
        EdgeData(vertex_data[5], vertex_data[0]),
        EdgeData(vertex_data[5], vertex_data[1]),
        EdgeData(vertex_data[6], vertex_data[0]),
        EdgeData(vertex_data[6], vertex_data[1]),

        EdgeData(vertex_data[7], vertex_data[1]),
        EdgeData(vertex_data[7], vertex_data[2]),
        EdgeData(vertex_data[8], vertex_data[1]),
        EdgeData(vertex_data[8], vertex_data[2]),

        EdgeData(vertex_data[9], vertex_data[2]),
        EdgeData(vertex_data[9], vertex_data[3]),
        EdgeData(vertex_data[10], vertex_data[2]),
        EdgeData(vertex_data[10], vertex_data[3]),

        EdgeData(vertex_data[11], vertex_data[3]),
        EdgeData(vertex_data[11], vertex_data[4]),
        EdgeData(vertex_data[12], vertex_data[3]),
        EdgeData(vertex_data[12], vertex_data[4]),

        EdgeData(vertex_data[13], vertex_data[4]),
        EdgeData(vertex_data[13], vertex_data[0]),
        EdgeData(vertex_data[14], vertex_data[4]),
        EdgeData(vertex_data[14], vertex_data[0])

    )

    num_windows = 7
    num_window_types = 2


class SixPlusTwelve(Cage):
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
        _CageVertexData(-1, -1, 0),
        _CageVertexData(-1, 1, 0),
        _CageVertexData(1, -1, 0),
        _CageVertexData(1, 1, 0),
        _CageVertexData(0, 0, 1),
        _CageVertexData(0, 0, -1)
    )

    vertex_data = (
        *_vertex_data,
        _CageVertexData.init_at_center(
            _vertex_data[0], _vertex_data[1]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[1], _vertex_data[3]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[3], _vertex_data[2]
        ),

        _CageVertexData.init_at_center(
            _vertex_data[0], _vertex_data[2]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[4], _vertex_data[0]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[4], _vertex_data[1]
        ),

        _CageVertexData.init_at_center(
            _vertex_data[4], _vertex_data[2]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[4], _vertex_data[3]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[5], _vertex_data[0]
        ),

        _CageVertexData.init_at_center(
            _vertex_data[5], _vertex_data[1]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[5], _vertex_data[2]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[5], _vertex_data[3]
        )
    )

    edge_data = (
        EdgeData(vertex_data[6], vertex_data[0]),
        EdgeData(vertex_data[6], vertex_data[1]),
        EdgeData(vertex_data[7], vertex_data[1]),
        EdgeData(vertex_data[7], vertex_data[3]),
        EdgeData(vertex_data[8], vertex_data[3]),
        EdgeData(vertex_data[8], vertex_data[2]),

        EdgeData(vertex_data[9], vertex_data[0]),
        EdgeData(vertex_data[9], vertex_data[2]),
        EdgeData(vertex_data[10], vertex_data[4]),
        EdgeData(vertex_data[10], vertex_data[0]),
        EdgeData(vertex_data[11], vertex_data[4]),
        EdgeData(vertex_data[11], vertex_data[1]),

        EdgeData(vertex_data[12], vertex_data[4]),
        EdgeData(vertex_data[12], vertex_data[2]),
        EdgeData(vertex_data[13], vertex_data[4]),
        EdgeData(vertex_data[13], vertex_data[3]),
        EdgeData(vertex_data[14], vertex_data[5]),
        EdgeData(vertex_data[14], vertex_data[0]),

        EdgeData(vertex_data[15], vertex_data[5]),
        EdgeData(vertex_data[15], vertex_data[1]),
        EdgeData(vertex_data[16], vertex_data[5]),
        EdgeData(vertex_data[16], vertex_data[2]),
        EdgeData(vertex_data[17], vertex_data[5]),
        EdgeData(vertex_data[17], vertex_data[3])
    )

    num_windows = 8
    num_window_types = 1


class EightPlusSixteen(Cage):
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

    _x = 2
    _vertex_data = (
        _CageVertexData(-0.5*_x, 0.5*_x, -0.35*_x),
        _CageVertexData(-0.5*_x, -0.5*_x, -0.35*_x),
        _CageVertexData(0.5*_x, -0.5*_x, -0.35*_x),
        _CageVertexData(0.5*_x, 0.5*_x, -0.35*_x),

        _CageVertexData(-_x*np.sqrt(2)/2, 0, _x*0.35),
        _CageVertexData(0, -_x*np.sqrt(2)/2, _x*0.35),
        _CageVertexData(_x*np.sqrt(2)/2, 0, _x*0.35),
        _CageVertexData(0, _x*np.sqrt(2)/2, _x*0.35)
    )

    vertex_data = (
        *_vertex_data,

        _CageVertexData.init_at_center(
            _vertex_data[1], _vertex_data[5]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[2], _vertex_data[5]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[0], _vertex_data[4]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[1], _vertex_data[4]
        ),

        _CageVertexData.init_at_center(
            _vertex_data[2], _vertex_data[6]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[3], _vertex_data[6]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[0], _vertex_data[7]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[3], _vertex_data[7]
        ),

        _CageVertexData.init_at_center(
            _vertex_data[0], _vertex_data[1]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[1], _vertex_data[2]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[2], _vertex_data[3]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[3], _vertex_data[0]
        ),

        _CageVertexData.init_at_center(
            _vertex_data[4], _vertex_data[5]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[5], _vertex_data[6]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[6], _vertex_data[7]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[7], _vertex_data[4]
        ),

    )

    edge_data = (
        EdgeData(vertex_data[8], vertex_data[1]),
        EdgeData(vertex_data[8], vertex_data[5]),

        EdgeData(vertex_data[9], vertex_data[2]),
        EdgeData(vertex_data[9], vertex_data[5]),

        EdgeData(vertex_data[10], vertex_data[0]),
        EdgeData(vertex_data[10], vertex_data[4]),

        EdgeData(vertex_data[11], vertex_data[1]),
        EdgeData(vertex_data[11], vertex_data[4]),

        EdgeData(vertex_data[12], vertex_data[2]),
        EdgeData(vertex_data[12], vertex_data[6]),

        EdgeData(vertex_data[13], vertex_data[3]),
        EdgeData(vertex_data[13], vertex_data[6]),

        EdgeData(vertex_data[14], vertex_data[0]),
        EdgeData(vertex_data[14], vertex_data[7]),

        EdgeData(vertex_data[15], vertex_data[3]),
        EdgeData(vertex_data[15], vertex_data[7]),

        EdgeData(vertex_data[16], vertex_data[0]),
        EdgeData(vertex_data[16], vertex_data[1]),

        EdgeData(vertex_data[17], vertex_data[1]),
        EdgeData(vertex_data[17], vertex_data[2]),

        EdgeData(vertex_data[18], vertex_data[2]),
        EdgeData(vertex_data[18], vertex_data[3]),

        EdgeData(vertex_data[19], vertex_data[3]),
        EdgeData(vertex_data[19], vertex_data[0]),

        EdgeData(vertex_data[20], vertex_data[4]),
        EdgeData(vertex_data[20], vertex_data[5]),

        EdgeData(vertex_data[21], vertex_data[5]),
        EdgeData(vertex_data[21], vertex_data[6]),

        EdgeData(vertex_data[22], vertex_data[6]),
        EdgeData(vertex_data[22], vertex_data[7]),

        EdgeData(vertex_data[23], vertex_data[7]),
        EdgeData(vertex_data[23], vertex_data[4])
    )

    num_windows = 10
    num_window_types = 2


class TenPlusTwenty(Cage):
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

    _x = 2
    _vertex_data = (
        _CageVertexData(-_x, _x, -_x),
        _CageVertexData(-_x, -_x, -_x),
        _CageVertexData(_x, _x, -_x),
        _CageVertexData(_x, -_x, -_x),

        _CageVertexData(-_x, _x, _x),
        _CageVertexData(-_x, -_x, _x),
        _CageVertexData(_x, _x, _x),
        _CageVertexData(_x, -_x, _x),

        _CageVertexData(0, 0, _x*1.5),
        _CageVertexData(0, 0, -_x*1.5)
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
        ),

        _CageVertexData.init_at_center(
            _vertex_data[8], _vertex_data[4]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[8], _vertex_data[5]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[8], _vertex_data[6]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[8], _vertex_data[7]
        ),

        _CageVertexData.init_at_center(
            _vertex_data[9], _vertex_data[0]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[9], _vertex_data[1]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[9], _vertex_data[2]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[9], _vertex_data[3]
        )
    )

    edge_data = (
        EdgeData(vertex_data[10], vertex_data[0]),
        EdgeData(vertex_data[10], vertex_data[2]),
        EdgeData(vertex_data[11], vertex_data[0]),
        EdgeData(vertex_data[11], vertex_data[1]),

        EdgeData(vertex_data[12], vertex_data[1]),
        EdgeData(vertex_data[12], vertex_data[3]),
        EdgeData(vertex_data[13], vertex_data[2]),
        EdgeData(vertex_data[13], vertex_data[3]),

        EdgeData(vertex_data[14], vertex_data[4]),
        EdgeData(vertex_data[14], vertex_data[6]),
        EdgeData(vertex_data[15], vertex_data[4]),
        EdgeData(vertex_data[15], vertex_data[5]),

        EdgeData(vertex_data[16], vertex_data[5]),
        EdgeData(vertex_data[16], vertex_data[7]),
        EdgeData(vertex_data[17], vertex_data[6]),
        EdgeData(vertex_data[17], vertex_data[7]),

        EdgeData(vertex_data[18], vertex_data[0]),
        EdgeData(vertex_data[18], vertex_data[4]),
        EdgeData(vertex_data[19], vertex_data[1]),
        EdgeData(vertex_data[19], vertex_data[5]),

        EdgeData(vertex_data[20], vertex_data[2]),
        EdgeData(vertex_data[20], vertex_data[6]),
        EdgeData(vertex_data[21], vertex_data[3]),
        EdgeData(vertex_data[21], vertex_data[7]),

        EdgeData(vertex_data[22], vertex_data[8]),
        EdgeData(vertex_data[22], vertex_data[4]),
        EdgeData(vertex_data[23], vertex_data[8]),
        EdgeData(vertex_data[23], vertex_data[5]),

        EdgeData(vertex_data[24], vertex_data[8]),
        EdgeData(vertex_data[24], vertex_data[6]),
        EdgeData(vertex_data[25], vertex_data[8]),
        EdgeData(vertex_data[25], vertex_data[7]),

        EdgeData(vertex_data[26], vertex_data[9]),
        EdgeData(vertex_data[26], vertex_data[0]),
        EdgeData(vertex_data[27], vertex_data[9]),
        EdgeData(vertex_data[27], vertex_data[1]),

        EdgeData(vertex_data[28], vertex_data[9]),
        EdgeData(vertex_data[28], vertex_data[2]),
        EdgeData(vertex_data[29], vertex_data[9]),
        EdgeData(vertex_data[29], vertex_data[3])
    )

    num_windows = 12
    num_window_types = 2
