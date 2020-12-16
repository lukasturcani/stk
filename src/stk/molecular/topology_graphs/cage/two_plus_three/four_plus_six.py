"""
Four Plus Six
=============

"""

import numpy as np

from ..cage import Cage
from ..vertices import _LinearCageVertex, _NonLinearCageVertex
from ...topology_graph import Edge


class FourPlusSix(Cage):
    """
    Represents a cage topology graph.

    Nonlinear building blocks with three functional groups are
    required for this topology.

    Linear building blocks with two functional groups are required for
    this topology.

    When using a :class:`dict` for the `building_blocks` parameter,
    as in :ref:`cage-topology-graph-examples`:
    *Multi-Building Block Cage Construction*, a
    :class:`.BuildingBlock`, with the following number of functional
    groups, needs to be assigned to each of the following vertex ids:

        | 3-functional groups: 0 to 3
        | 2-functional groups: 4 to 9

    See :class:`.Cage` for more details and examples.

    """

    # Vertices of a tetrahdron so that origin is at the origin. Source:
    # http://tinyurl.com/lc262h8.
    _v0, _v1, _v2, _v3 = _vertex_prototypes = (
        _NonLinearCageVertex(0, [0, 0, np.sqrt(6)/2]),
        _NonLinearCageVertex(1, [-1, -np.sqrt(3)/3, -np.sqrt(6)/6]),
        _NonLinearCageVertex(2, [1, -np.sqrt(3)/3, -np.sqrt(6)/6]),
        _NonLinearCageVertex(3, [0, 2*np.sqrt(3)/3, -np.sqrt(6)/6]),
    )

    _vertex_prototypes = (
        *_vertex_prototypes,
        _LinearCageVertex.init_at_center(4, (_v0, _v1)),
        _LinearCageVertex.init_at_center(5, (_v0, _v2)),
        _LinearCageVertex.init_at_center(6, (_v0, _v3)),
        _LinearCageVertex.init_at_center(7, (_v1, _v2)),
        _LinearCageVertex.init_at_center(8, (_v1, _v3)),
        _LinearCageVertex.init_at_center(9, (_v2, _v3)),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[4], _vertex_prototypes[0]),
        Edge(1, _vertex_prototypes[4], _vertex_prototypes[1]),
        Edge(2, _vertex_prototypes[5], _vertex_prototypes[0]),
        Edge(3, _vertex_prototypes[5], _vertex_prototypes[2]),
        Edge(4, _vertex_prototypes[6], _vertex_prototypes[0]),
        Edge(5, _vertex_prototypes[6], _vertex_prototypes[3]),
        Edge(6, _vertex_prototypes[7], _vertex_prototypes[1]),
        Edge(7, _vertex_prototypes[7], _vertex_prototypes[2]),
        Edge(8, _vertex_prototypes[8], _vertex_prototypes[1]),
        Edge(9, _vertex_prototypes[8], _vertex_prototypes[3]),
        Edge(10, _vertex_prototypes[9], _vertex_prototypes[2]),
        Edge(11, _vertex_prototypes[9], _vertex_prototypes[3])
    )

    _num_windows = 4
    _num_window_types = 1
