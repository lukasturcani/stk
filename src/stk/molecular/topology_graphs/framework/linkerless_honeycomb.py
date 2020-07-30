"""
Linkerless Honeycomb
====================

"""

import numpy as np

from .framework import Framework
from .vertices import _NonLinearFrameworkVertex
from ..topology_graph import Edge


class LinkerlessHoneycomb(Framework):
    """
    Represents a honeycomb framework topology graph.

    Building blocks with three functional groups are required
    for this topology graph.

    See :class:`.Framework` for more details and examples.

    """

    _lattice_constants = _a, _b, _c = (
        np.array([1., 0., 0.]),
        np.array([0.5, 0.866, 0.]),
        np.array([0., 0., 5/1.7321]),
    )

    _vertex_prototypes = (
        _NonLinearFrameworkVertex(0, (1/3)*_a + (1/3)*_b + (1/2)*_c),
        _NonLinearFrameworkVertex(1, (2/3)*_a + (2/3)*_b + (1/2)*_c),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[1]),
        Edge(
            id=1,
            vertex1=_vertex_prototypes[0],
            vertex2=_vertex_prototypes[1],
            periodicity=(-1, 0, 0),
        ),
        Edge(
            id=2,
            vertex1=_vertex_prototypes[0],
            vertex2=_vertex_prototypes[1],
            periodicity=(0, -1, 0),
        ),
    )
