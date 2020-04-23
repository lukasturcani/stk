"""
Hexagonal
=========

"""

import numpy as np


from .cof import Cof
from .vertices import _LinearCofVertex, _NonLinearCofVertex
from ..topology_graph import Edge


class Hexagonal(Cof):
    """
    Represents a hexagonal COF topology graph.

    Building blocks with six and two functional groups are required
    for this topology graph.

    See :class:`.Cof` for more details and examples.

    """

    _lattice_constants = _a, _b, _c = (
        np.array([1., 0., 0.]),
        np.array([0.5, 0.866, 0]),
        np.array([0, 0, 5/1.7321])
    )

    _vertex_prototypes = (
        _NonLinearCofVertex(0, (1/4)*_a + (1/4)*_b + (1/2)*_c),
        _NonLinearCofVertex(1, (1/4)*_a + (3/4)*_b + (1/2)*_c),
        _NonLinearCofVertex(2, (3/4)*_a + (1/4)*_b + (1/2)*_c),
        _NonLinearCofVertex(3, (3/4)*_a + (3/4)*_b + (1/2)*_c),
    )

    _vertex_prototypes = (
        *_vertex_prototypes,
        _LinearCofVertex.init_at_center(
            id=4,
            vertices=(_vertex_prototypes[0], _vertex_prototypes[1]),
        ),
        _LinearCofVertex.init_at_center(
            id=5,
            vertices=(_vertex_prototypes[0], _vertex_prototypes[2]),
        ),
        _LinearCofVertex.init_at_center(
            id=6,
            vertices=(_vertex_prototypes[1], _vertex_prototypes[2]),
        ),
        _LinearCofVertex.init_at_center(
            id=7,
            vertices=(_vertex_prototypes[1], _vertex_prototypes[3]),
        ),
        _LinearCofVertex.init_at_center(
            id=8,
            vertices=(_vertex_prototypes[2], _vertex_prototypes[3]),
        ),
        _LinearCofVertex.init_at_shifted_center(
            id=9,
            vertices=(_vertex_prototypes[0], _vertex_prototypes[2]),
            cell_shifts=((0, 0, 0), (-1, 0, 0)),
            lattice_constants=_lattice_constants,
        ),
        _LinearCofVertex.init_at_shifted_center(
            id=10,
            vertices=(_vertex_prototypes[0], _vertex_prototypes[1]),
            cell_shifts=((0, 0, 0), (0, -1, 0)),
            lattice_constants=_lattice_constants,
        ),
        _LinearCofVertex.init_at_shifted_center(
            id=11,
            vertices=(_vertex_prototypes[0], _vertex_prototypes[3]),
            cell_shifts=((0, 0, 0), (0, -1, 0)),
            lattice_constants=_lattice_constants,
        ),
        _LinearCofVertex.init_at_shifted_center(
            id=12,
            vertices=(_vertex_prototypes[2], _vertex_prototypes[1]),
            cell_shifts=((0, 0, 0), (1, -1, 0)),
            lattice_constants=_lattice_constants,
        ),
        _LinearCofVertex.init_at_shifted_center(
            id=13,
            vertices=(_vertex_prototypes[2], _vertex_prototypes[3]),
            cell_shifts=((0, 0, 0), (0, -1, 0)),
            lattice_constants=_lattice_constants,
        ),
        _LinearCofVertex.init_at_shifted_center(
            id=14,
            vertices=(_vertex_prototypes[1], _vertex_prototypes[3]),
            cell_shifts=((0, 0, 0), (-1, 0, 0)),
            lattice_constants=_lattice_constants,
        ),
        _LinearCofVertex.init_at_shifted_center(
            id=15,
            vertices=(_vertex_prototypes[3], _vertex_prototypes[0]),
            cell_shifts=((0, 0, 0), (1, 0, 0)),
            lattice_constants=_lattice_constants,
        )
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[4], _vertex_prototypes[0]),
        Edge(1, _vertex_prototypes[4], _vertex_prototypes[1]),

        Edge(2, _vertex_prototypes[5], _vertex_prototypes[0]),
        Edge(3, _vertex_prototypes[5], _vertex_prototypes[2]),

        Edge(4, _vertex_prototypes[6], _vertex_prototypes[1]),
        Edge(5, _vertex_prototypes[6], _vertex_prototypes[2]),

        Edge(6, _vertex_prototypes[7], _vertex_prototypes[1]),
        Edge(7, _vertex_prototypes[7], _vertex_prototypes[3]),

        Edge(8, _vertex_prototypes[8], _vertex_prototypes[2]),
        Edge(9, _vertex_prototypes[8], _vertex_prototypes[3]),

        Edge(10, _vertex_prototypes[9], _vertex_prototypes[0]),
        Edge(
            id=11,
            vertex1=_vertex_prototypes[9],
            vertex2=_vertex_prototypes[2],
            periodicity=(-1, 0, 0),
        ),

        Edge(12, _vertex_prototypes[10], _vertex_prototypes[0]),
        Edge(
            id=13,
            vertex1=_vertex_prototypes[10],
            vertex2=_vertex_prototypes[1],
            periodicity=(0, -1, 0),
        ),

        Edge(14, _vertex_prototypes[11], _vertex_prototypes[0]),
        Edge(
            id=15,
            vertex1=_vertex_prototypes[11],
            vertex2=_vertex_prototypes[3],
            periodicity=(0, -1, 0),
        ),

        Edge(16, _vertex_prototypes[12], _vertex_prototypes[2]),
        Edge(
            id=17,
            vertex1=_vertex_prototypes[12],
            vertex2=_vertex_prototypes[1],
            periodicity=(1, -1, 0),
        ),

        Edge(18, _vertex_prototypes[13], _vertex_prototypes[2]),
        Edge(
            id=19,
            vertex1=_vertex_prototypes[13],
            vertex2=_vertex_prototypes[3],
            periodicity=(0, -1, 0),
        ),

        Edge(20, _vertex_prototypes[14], _vertex_prototypes[1]),
        Edge(
            id=21,
            vertex1=_vertex_prototypes[14],
            vertex2=_vertex_prototypes[3],
            periodicity=(-1, 0, 0),
        ),

        Edge(22, _vertex_prototypes[15], _vertex_prototypes[3]),
        Edge(
            id=23,
            vertex1=_vertex_prototypes[15],
            vertex2=_vertex_prototypes[0],
            periodicity=(1, 0, 0),
        ),
    )
