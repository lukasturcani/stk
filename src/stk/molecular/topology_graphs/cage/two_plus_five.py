"""
Defines cages made from building blocks with 2 and 5 functional groups.

"""

from scipy.constants import golden

from .base import Cage, _CageVertexData
from ..topology_graph import EdgeData


class TwelvePlusThirty(Cage):
    """
    Represents a icosahedron cage topology graph.

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

    # Vertices of a regular origin-centred icosahedron
    # Source: http://eusebeia.dyndns.org/4d/icosahedron
    _vertex_data = (
        _CageVertexData(0, 1, golden),
        _CageVertexData(0, -1, golden),
        _CageVertexData(0, 1, -golden),
        _CageVertexData(0, -1, -golden),
        _CageVertexData(1, golden, 0),
        _CageVertexData(-1, golden, 0),
        _CageVertexData(1, -golden, 0),
        _CageVertexData(-1, -golden, 0),
        _CageVertexData(golden, 0, 1),
        _CageVertexData(-golden, 0, 1),
        _CageVertexData(golden, 0, -1),
        _CageVertexData(-golden, 0, -1)
    )

    vertex_data = (
        *_vertex_data,
        _CageVertexData.init_at_center(
            _vertex_data[0], _vertex_data[1]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[0], _vertex_data[9]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[0], _vertex_data[5]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[0], _vertex_data[4]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[0], _vertex_data[8]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[8], _vertex_data[1]),
        _CageVertexData.init_at_center(
            _vertex_data[1], _vertex_data[9]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[9], _vertex_data[5]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[5], _vertex_data[4]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[4], _vertex_data[8]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[5], _vertex_data[2]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[5], _vertex_data[11]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[9], _vertex_data[11]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[9], _vertex_data[7]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[1], _vertex_data[7]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[1], _vertex_data[6]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[8], _vertex_data[6]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[8], _vertex_data[10]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[4], _vertex_data[10]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[4], _vertex_data[2]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[2], _vertex_data[11]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[11], _vertex_data[7]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[7], _vertex_data[6]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[6], _vertex_data[10]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[10], _vertex_data[2]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[2], _vertex_data[3]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[11], _vertex_data[3]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[7], _vertex_data[3]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[6], _vertex_data[3]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[10], _vertex_data[3]
        )
    )

    edge_data = (
        EdgeData(vertex_data[12], vertex_data[0]),
        EdgeData(vertex_data[12], vertex_data[1]),
        EdgeData(vertex_data[13], vertex_data[0]),
        EdgeData(vertex_data[13], vertex_data[9]),

        EdgeData(vertex_data[14], vertex_data[0]),
        EdgeData(vertex_data[14], vertex_data[5]),
        EdgeData(vertex_data[15], vertex_data[0]),
        EdgeData(vertex_data[15], vertex_data[4]),

        EdgeData(vertex_data[16], vertex_data[0]),
        EdgeData(vertex_data[16], vertex_data[8]),
        EdgeData(vertex_data[17], vertex_data[8]),
        EdgeData(vertex_data[17], vertex_data[1]),

        EdgeData(vertex_data[18], vertex_data[1]),
        EdgeData(vertex_data[18], vertex_data[9]),
        EdgeData(vertex_data[19], vertex_data[9]),
        EdgeData(vertex_data[19], vertex_data[5]),

        EdgeData(vertex_data[20], vertex_data[5]),
        EdgeData(vertex_data[20], vertex_data[4]),
        EdgeData(vertex_data[21], vertex_data[4]),
        EdgeData(vertex_data[21], vertex_data[8]),

        EdgeData(vertex_data[22], vertex_data[5]),
        EdgeData(vertex_data[22], vertex_data[2]),
        EdgeData(vertex_data[23], vertex_data[5]),
        EdgeData(vertex_data[23], vertex_data[11]),

        EdgeData(vertex_data[24], vertex_data[9]),
        EdgeData(vertex_data[24], vertex_data[11]),
        EdgeData(vertex_data[25], vertex_data[9]),
        EdgeData(vertex_data[25], vertex_data[7]),

        EdgeData(vertex_data[26], vertex_data[1]),
        EdgeData(vertex_data[26], vertex_data[7]),
        EdgeData(vertex_data[27], vertex_data[1]),
        EdgeData(vertex_data[27], vertex_data[6]),

        EdgeData(vertex_data[28], vertex_data[8]),
        EdgeData(vertex_data[28], vertex_data[6]),
        EdgeData(vertex_data[29], vertex_data[8]),
        EdgeData(vertex_data[29], vertex_data[10]),

        EdgeData(vertex_data[30], vertex_data[4]),
        EdgeData(vertex_data[30], vertex_data[10]),
        EdgeData(vertex_data[31], vertex_data[4]),
        EdgeData(vertex_data[31], vertex_data[2]),

        EdgeData(vertex_data[32], vertex_data[2]),
        EdgeData(vertex_data[32], vertex_data[11]),
        EdgeData(vertex_data[33], vertex_data[11]),
        EdgeData(vertex_data[33], vertex_data[7]),

        EdgeData(vertex_data[34], vertex_data[7]),
        EdgeData(vertex_data[34], vertex_data[6]),
        EdgeData(vertex_data[35], vertex_data[6]),
        EdgeData(vertex_data[35], vertex_data[10]),

        EdgeData(vertex_data[36], vertex_data[10]),
        EdgeData(vertex_data[36], vertex_data[2]),
        EdgeData(vertex_data[37], vertex_data[2]),
        EdgeData(vertex_data[37], vertex_data[3]),

        EdgeData(vertex_data[38], vertex_data[11]),
        EdgeData(vertex_data[38], vertex_data[3]),
        EdgeData(vertex_data[39], vertex_data[7]),
        EdgeData(vertex_data[39], vertex_data[3]),

        EdgeData(vertex_data[40], vertex_data[6]),
        EdgeData(vertex_data[40], vertex_data[3]),
        EdgeData(vertex_data[41], vertex_data[10]),
        EdgeData(vertex_data[41], vertex_data[3])
    )

    num_windows = 20
    num_window_types = 1
