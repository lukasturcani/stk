"""
Linkerless Honeycomb
====================

"""

import numpy as np

from .cof import Cof
from .vertices import _NonLinearCofVertex
from ..topology_graph import Edge


class LinkerlessHoneycomb(Cof):
    """
    Represents a honeycomb COF topology graph.

    Building blocks with three functional groups are required
    for this topology graph.

    When using a :class:`dict` for the `building_blocks` parameter,
    as in :ref:`cof-topology-graph-examples`:
    *Multi-Building Block COF Construction*, a
    :class:`.BuildingBlock`, with the following number of functional
    groups, needs to be assigned to each of the following vertex ids:

        | 3-functional groups: 0 to 1

    See :class:`.Cof` for more details and examples.

    """

    _lattice_constants = _a, _b, _c = (
        np.array([1., 0., 0.]),
        np.array([0.5, 0.866, 0.]),
        np.array([0., 0., 5/1.7321]),
    )

    _vertex_prototypes = (
        _NonLinearCofVertex(0, (1/3)*_a + (1/3)*_b + (1/2)*_c),
        _NonLinearCofVertex(1, (2/3)*_a + (2/3)*_b + (1/2)*_c),
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
