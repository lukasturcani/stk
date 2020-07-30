"""
Square
======

"""

import numpy as np


from .framework import Framework
from .vertices import _LinearFrameworkVertex, _NonLinearFrameworkVertex
from ..topology_graph import Edge


class Square(Framework):
    """
    Represents a sqaure framework topology graph.

    Building blocks with four and two functional groups are required
    for this topology graph.

    See :class:`.Framework` for more details and examples.

    """

    _lattice_constants = _a, _b, _c = (
        np.array([1., 0., 0.]),
        np.array([0., 1., 0.]),
        np.array([0., 0., 1.])
    )

    _vertex_prototypes = (
        _NonLinearFrameworkVertex(0, (0.5)*_a + (0.5)*_b + (0.5)*_c),
    )
    _vertex_prototypes = (
        *_vertex_prototypes,
        _LinearFrameworkVertex.init_at_shifted_center(
            id=1,
            vertices=(_vertex_prototypes[0], _vertex_prototypes[0]),
            cell_shifts=((0, 0, 0), (1, 0, 0)),
            lattice_constants=_lattice_constants,
        ),
        _LinearFrameworkVertex.init_at_shifted_center(
            id=2,
            vertices=(_vertex_prototypes[0], _vertex_prototypes[0]),
            cell_shifts=((0, 0, 0), (0, 1, 0)),
            lattice_constants=_lattice_constants,
        ),

    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[1], _vertex_prototypes[0]),
        Edge(
            id=1,
            vertex1=_vertex_prototypes[1],
            vertex2=_vertex_prototypes[0],
            periodicity=(1, 0, 0),
        ),
        Edge(2, _vertex_prototypes[2], _vertex_prototypes[0]),
        Edge(
            id=3,
            vertex1=_vertex_prototypes[2],
            vertex2=_vertex_prototypes[0],
            periodicity=(0, 1, 0),
        ),
    )
