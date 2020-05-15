"""
M3L6
====

"""

import numpy as np

from ..cage import Cage
from ..vertices import _NonLinearCageVertex, _LinearCageVertex
from ...topology_graph import Edge


class M3L6(Cage):
    """
    Represents a cage topology graph.

    Metal building blocks with four functional groups are
    required for this topology.

    Ligand building blocks with two functional groups are required for
    this topology.

    When using a :class:`dict` for initialization, a
    :class:`.BuildingBlock` needs to be assigned to each of the
    following numbers:

        | metals: (0, 1, 2)
        | ligands: (3, 4, 5, 6, 7, 8)

    See :class:`.Cage` for more details and examples.

    """

    _R, _theta = 1, 0

    _vertex_prototypes = (
        _NonLinearCageVertex(
            id=0,
            position=[_R*np.cos(_theta), _R*np.sin(_theta), 0]
        ),
        _NonLinearCageVertex(
            id=1,
            position=[
                _R*np.cos(_theta+(4*np.pi/3)),
                _R*np.sin(_theta+(4*np.pi/3)),
                0
            ]
        ),
        _NonLinearCageVertex(
            id=2,
            position=[
                _R*np.cos(_theta+(2*np.pi/3)),
                _R*np.sin(_theta+(2*np.pi/3)),
                0
            ]
        ),

        _LinearCageVertex(
            id=3,
            position=[
                _R*np.cos((_theta+np.pi/4)),
                _R*np.sin((_theta+np.pi/4)),
                1
            ],
            use_neighbor_placement=False
        ),
        _LinearCageVertex(
            id=4,
            position=[
                _R*np.cos((_theta+1*np.pi/3)),
                _R*np.sin((_theta+1*np.pi/3)),
                -1
            ],
            use_neighbor_placement=False
        ),

        _LinearCageVertex(
            id=5,
            position=[
                _R*np.cos((_theta+1*np.pi/3)+(4*np.pi/3)),
                _R*np.sin((_theta+1*np.pi/3)+(4*np.pi/3)),
                1
            ],
            use_neighbor_placement=False
        ),
        _LinearCageVertex(
            id=6,
            position=[
                _R*np.cos((_theta+1*np.pi/3)+(4*np.pi/3)),
                _R*np.sin((_theta+1*np.pi/3)+(4*np.pi/3)),
                -1
            ],
            use_neighbor_placement=False
        ),

        _LinearCageVertex(
            id=7,
            position=[
                _R*np.cos((_theta+1*np.pi/3)+(2*np.pi/3)),
                _R*np.sin((_theta+1*np.pi/3)+(2*np.pi/3)),
                1
            ],
            use_neighbor_placement=False
        ),
        _LinearCageVertex(
            id=8,
            position=[
                _R*np.cos((_theta+1*np.pi/3)+(2*np.pi/3)),
                _R*np.sin((_theta+1*np.pi/3)+(2*np.pi/3)),
                -1
            ],
            use_neighbor_placement=False
        ),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[3]),
        Edge(1, _vertex_prototypes[0], _vertex_prototypes[4]),
        Edge(2, _vertex_prototypes[0], _vertex_prototypes[5]),
        Edge(3, _vertex_prototypes[0], _vertex_prototypes[6]),

        Edge(4, _vertex_prototypes[1], _vertex_prototypes[5]),
        Edge(5, _vertex_prototypes[1], _vertex_prototypes[6]),
        Edge(6, _vertex_prototypes[1], _vertex_prototypes[7]),
        Edge(7, _vertex_prototypes[1], _vertex_prototypes[8]),

        Edge(8, _vertex_prototypes[2], _vertex_prototypes[3]),
        Edge(9, _vertex_prototypes[2], _vertex_prototypes[4]),
        Edge(10, _vertex_prototypes[2], _vertex_prototypes[7]),
        Edge(11, _vertex_prototypes[2], _vertex_prototypes[8]),
    )

    _num_windows = 2
    _num_window_types = 1
