"""
Kagome
======

"""

import numpy as np

from .cof import Cof
from .vertices import LinearCofVertex, NonLinearCofVertex
from ..topology_graph import Edge


class Kagome(Cof):
    """
    Represents a kagome COF topology graph.

    Building blocks with four and two functional groups are required
    for this topology graph.

    When using a :class:`dict` for the `building_blocks` parameter,
    as in :ref:`cof-topology-graph-examples`:
    *Multi-Building Block COF Construction*, a
    :class:`.BuildingBlock`, with the following number of functional
    groups, needs to be assigned to each of the following vertex ids:

        | 4-functional groups: 0 to 2
        | 2-functional groups: 3 to 8

    See :class:`.Cof` for more details and examples.

    """

    _lattice_constants = _a, _b, _c = (
        np.array([1., 0., 0.]),
        np.array([0.5, 0.866, 0.]),
        np.array([0., 0., 5/1.7321])
    )

    _vertex_prototypes = (
        NonLinearCofVertex(0, (1/4)*_a + (3/4)*_b + (0.5)*_c),
        NonLinearCofVertex(1, (3/4)*_a + (3/4)*_b + (1/2)*_c),
        NonLinearCofVertex(2, (3/4)*_a + (1/4)*_b + (1/2)*_c),
    )

    _vertex_prototypes = (
        *_vertex_prototypes,
        LinearCofVertex.init_at_center(
            id=3,
            vertices=(_vertex_prototypes[0], _vertex_prototypes[1]),
        ),
        LinearCofVertex.init_at_center(
            id=4,
            vertices=(_vertex_prototypes[0], _vertex_prototypes[2]),
        ),
        LinearCofVertex.init_at_center(
            id=5,
            vertices=(_vertex_prototypes[1], _vertex_prototypes[2]),
        ),
        LinearCofVertex.init_at_shifted_center(
            id=6,
            vertices=(_vertex_prototypes[0], _vertex_prototypes[1]),
            cell_shifts=((0, 0, 0), (-1, 0, 0)),
            lattice_constants=_lattice_constants
        ),
        LinearCofVertex.init_at_shifted_center(
            id=7,
            vertices=(_vertex_prototypes[0], _vertex_prototypes[2]),
            cell_shifts=((0, 0, 0), (-1, 1, 0)),
            lattice_constants=_lattice_constants
        ),
        LinearCofVertex.init_at_shifted_center(
            id=8,
            vertices=(_vertex_prototypes[1], _vertex_prototypes[2]),
            cell_shifts=((0, 0, 0), (0, 1, 0)),
            lattice_constants=_lattice_constants
        ),

    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[3], _vertex_prototypes[0]),
        Edge(1, _vertex_prototypes[3], _vertex_prototypes[1]),

        Edge(2, _vertex_prototypes[4], _vertex_prototypes[0]),
        Edge(3, _vertex_prototypes[4], _vertex_prototypes[2]),

        Edge(4, _vertex_prototypes[5], _vertex_prototypes[1]),
        Edge(5, _vertex_prototypes[5], _vertex_prototypes[2]),

        Edge(6, _vertex_prototypes[6], _vertex_prototypes[0]),
        Edge(
            id=7,
            vertex1=_vertex_prototypes[6],
            vertex2=_vertex_prototypes[1],
            periodicity=(-1, 0, 0),
        ),

        Edge(8, _vertex_prototypes[7], _vertex_prototypes[0]),
        Edge(
            id=9,
            vertex1=_vertex_prototypes[7],
            vertex2=_vertex_prototypes[2],
            periodicity=(-1, 1, 0),
        ),

        Edge(10, _vertex_prototypes[8], _vertex_prototypes[1]),
        Edge(
            id=11,
            vertex1=_vertex_prototypes[8],
            vertex2=_vertex_prototypes[2],
            periodicity=(0, 1, 0),
        ),
    )
