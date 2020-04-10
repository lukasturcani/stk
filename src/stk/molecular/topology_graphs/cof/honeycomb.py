"""
Honeycomb
=========

"""

import numpy as np


from .cof import Cof
from .vertices import _LinearCofVertex, _NonLinearCofVertex
from ..topology_graph import Edge


class Honeycomb(Cof):
    """
    Represents a honeycomb COF topology graph.

    Building blocks with three and two functional groups are required
    for this topology graph.

    See :class:`.Cof` for more details and examples.

    """

    _lattice_constants = _a, _b, _c = (
        np.array([1., 0., 0.]),
        np.array([0.5, 0.866, 0]),
        np.array([0, 0, 5/1.7321])
    )

    _vertex_prototypes = (
        _NonLinearCofVertex(0, (1/3)*_a + (1/3)*_b + (1/2)*_c),
        _NonLinearCofVertex(1, (2/3)*_a + (2/3)*_b + (1/2)*_c),
    )

    _vertex_prototypes = (
        *_vertex_prototypes,
        _LinearCofVertex.init_at_center(
            id=2,
            vertices=(_vertex_prototypes[0], _vertex_prototypes[1]),
        ),
        _LinearCofVertex.init_at_shifted_center(
            id=3,
            vertices=(_vertex_prototypes[0], _vertex_prototypes[1]),
            cell_shifts=((0, 0, 0), (0, -1, 0)),
            lattice_constants=_lattice_constants,
        ),
        _LinearCofVertex.init_at_shifted_center(
            id=4,
            vertices=(_vertex_prototypes[0], _vertex_prototypes[1]),
            cell_shifts=((0, 0, 0), (-1, 0, 0)),
            lattice_constants=_lattice_constants,
        )
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[2], _vertex_prototypes[0]),
        Edge(1, _vertex_prototypes[2], _vertex_prototypes[1]),

        Edge(2, _vertex_prototypes[3], _vertex_prototypes[0]),
        Edge(
            id=3,
            vertex1=_vertex_prototypes[3],
            vertex2=_vertex_prototypes[1],
            periodicity=(0, -1, 0),
        ),

        Edge(4, _vertex_prototypes[4], _vertex_prototypes[0]),
        Edge(
            id=5,
            vertex1=_vertex_prototypes[4],
            vertex2=_vertex_prototypes[1],
            periodicity=(-1, 0, 0),
        )
    )
