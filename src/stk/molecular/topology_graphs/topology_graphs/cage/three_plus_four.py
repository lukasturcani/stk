"""
Defines cages made from building blocks with 3 and 4 functional groups.

"""

from .base import Cage, _CageVertexData
from ..topology_graph import EdgeData


class SixPlusEight(Cage):
    """
    Represents a cage topology graph.

    Building blocks with three and four functional groups are required
    for this topology graph.

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
    _vertex_data = (
        _CageVertexData(-_x, _x, 0),
        _CageVertexData(-_x, -_x, 0),
        _CageVertexData(_x, _x, 0),
        _CageVertexData(_x, -_x, 0),

        _CageVertexData(0, 0, _x),
        _CageVertexData(0, 0, -_x),
    )

    vertex_data = (
        *_vertex_data,
        _CageVertexData.init_at_center(
            _vertex_data[0],
            _vertex_data[4],
            _vertex_data[1]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[1],
            _vertex_data[4],
            _vertex_data[3]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[4],
            _vertex_data[3],
            _vertex_data[2]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[4],
            _vertex_data[2],
            _vertex_data[0]
        ),

        _CageVertexData.init_at_center(
            _vertex_data[0],
            _vertex_data[5],
            _vertex_data[1]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[5],
            _vertex_data[1],
            _vertex_data[3]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[3],
            _vertex_data[5],
            _vertex_data[2]
        ),
        _CageVertexData.init_at_center(
            _vertex_data[2],
            _vertex_data[5],
            _vertex_data[0]
        )
    )

    edge_data = (
        EdgeData(vertex_data[6], vertex_data[0]),
        EdgeData(vertex_data[6], vertex_data[4]),
        EdgeData(vertex_data[6], vertex_data[1]),

        EdgeData(vertex_data[7], vertex_data[1]),
        EdgeData(vertex_data[7], vertex_data[4]),
        EdgeData(vertex_data[7], vertex_data[3]),

        EdgeData(vertex_data[8], vertex_data[4]),
        EdgeData(vertex_data[8], vertex_data[3]),
        EdgeData(vertex_data[8], vertex_data[2]),

        EdgeData(vertex_data[9], vertex_data[4]),
        EdgeData(vertex_data[9], vertex_data[2]),
        EdgeData(vertex_data[9], vertex_data[0]),

        EdgeData(vertex_data[10], vertex_data[0]),
        EdgeData(vertex_data[10], vertex_data[5]),
        EdgeData(vertex_data[10], vertex_data[1]),

        EdgeData(vertex_data[11], vertex_data[5]),
        EdgeData(vertex_data[11], vertex_data[1]),
        EdgeData(vertex_data[11], vertex_data[3]),

        EdgeData(vertex_data[12], vertex_data[3]),
        EdgeData(vertex_data[12], vertex_data[5]),
        EdgeData(vertex_data[12], vertex_data[2]),

        EdgeData(vertex_data[13], vertex_data[2]),
        EdgeData(vertex_data[13], vertex_data[5]),
        EdgeData(vertex_data[13], vertex_data[0]),

    )

    num_windows = 12
    num_window_types = 1
