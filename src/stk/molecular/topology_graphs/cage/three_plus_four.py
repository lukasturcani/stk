"""
Defines cages made from building blocks with 3 and 4 functional groups.

"""

from .base import Cage, _CageVertex
from ..topology_graph import Edge


class SixPlusEight(Cage):
    """
    Represents a cage topology graph.

    Building blocks with three and four functional groups are required
    for this topology graph.

    See :class:`.Cage` for more details and examples.

    Attributes
    ----------
    vertices : :class:`tuple` of :class:`.Vertex`
        The vertices which make up the topology graph.

    edges : :class:`tuple` of :class:`.Edge`
        The edges which make up the topology graph.

    """

    _x = 1
    _vertices = (
        _CageVertex(-_x, _x, 0),
        _CageVertex(-_x, -_x, 0),
        _CageVertex(_x, _x, 0),
        _CageVertex(_x, -_x, 0),

        _CageVertex(0, 0, _x),
        _CageVertex(0, 0, -_x),
    )

    vertices = (
        *_vertices,
        _CageVertex.init_at_center(
            _vertices[0],
            _vertices[4],
            _vertices[1]
        ),
        _CageVertex.init_at_center(
            _vertices[1],
            _vertices[4],
            _vertices[3]
        ),
        _CageVertex.init_at_center(
            _vertices[4],
            _vertices[3],
            _vertices[2]
        ),
        _CageVertex.init_at_center(
            _vertices[4],
            _vertices[2],
            _vertices[0]
        ),

        _CageVertex.init_at_center(
            _vertices[0],
            _vertices[5],
            _vertices[1]
        ),
        _CageVertex.init_at_center(
            _vertices[5],
            _vertices[1],
            _vertices[3]
        ),
        _CageVertex.init_at_center(
            _vertices[3],
            _vertices[5],
            _vertices[2]
        ),
        _CageVertex.init_at_center(
            _vertices[2],
            _vertices[5],
            _vertices[0]
        )
    )

    edges = (
        Edge(vertices[6], vertices[0]),
        Edge(vertices[6], vertices[4]),
        Edge(vertices[6], vertices[1]),

        Edge(vertices[7], vertices[1]),
        Edge(vertices[7], vertices[4]),
        Edge(vertices[7], vertices[3]),

        Edge(vertices[8], vertices[4]),
        Edge(vertices[8], vertices[3]),
        Edge(vertices[8], vertices[2]),

        Edge(vertices[9], vertices[4]),
        Edge(vertices[9], vertices[2]),
        Edge(vertices[9], vertices[0]),

        Edge(vertices[10], vertices[0]),
        Edge(vertices[10], vertices[5]),
        Edge(vertices[10], vertices[1]),

        Edge(vertices[11], vertices[5]),
        Edge(vertices[11], vertices[1]),
        Edge(vertices[11], vertices[3]),

        Edge(vertices[12], vertices[3]),
        Edge(vertices[12], vertices[5]),
        Edge(vertices[12], vertices[2]),

        Edge(vertices[13], vertices[2]),
        Edge(vertices[13], vertices[5]),
        Edge(vertices[13], vertices[0]),

    )

    num_windows = 12
    num_window_types = 1
