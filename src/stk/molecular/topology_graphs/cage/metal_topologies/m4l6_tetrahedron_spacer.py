"""
M4L6 Tetrahedron with Spacer
============================

"""

import numpy as np

from ..cage import Cage
from ..vertices import _NonLinearCageVertex, _LinearCageVertex
from ...topology_graph import Edge


class M4L6TetrahedronSpacer(Cage):
    """
    Represents a cage topology graph.

    This topology places a ditopic spacer between the vertices of the
    tetrahderon.

    Metal-based building blocks with three functional groups are
    required for this topology.

    Linker building blocks with two functional groups are required for
    this topology.

    When using a :class:`dict` for initialization, a
    :class:`.BuildingBlock` needs to be assigned to each of the
    following numbers:

        | metals: (0, 1, 2, 3)
        | linkers: (4, 5, 6, 7, 8, 9)

    Examples
    --------
    *Building Metal-Organic Tetrahedron*

    Many metal-organic cages are built using a process called
    subcomponent self-assembly, which is a complex chemical process
    that occurs in solution. In :class:`.Cage`, an example is
    provided of an alchemical approach using
    :class:`.M4L6TetrahedronSpacer` to
    build these types of cages. It is alchemical because the bonds
    formed during construction are not the same as the experimental
    reaction. Instead of forming bonds at the metal centre, we create
    bonds between disconnected ligands.

    The :class:`.M4L6TetrahedronSpacer` topology is provided for cases
    where the linkers cannot be disconnected symmetrically. However,
    in the case that the linker can be disconnected
    in a symmetrical fashion, the
    :class:`.M4L6Tetrahedron` topology can be used.

    """

    _vertex_prototypes = (
        _NonLinearCageVertex(0, [0, 0, np.sqrt(6)/2]),
        _NonLinearCageVertex(1, [-1, -np.sqrt(3)/3, -np.sqrt(6)/6]),
        _NonLinearCageVertex(2, [1, -np.sqrt(3)/3, -np.sqrt(6)/6]),
        _NonLinearCageVertex(3, [0, 2*np.sqrt(3)/3, -np.sqrt(6)/6]),
    )

    _vertex_prototypes = (
        *_vertex_prototypes,

        _LinearCageVertex.init_at_center(
            id=4,
            vertices=(_vertex_prototypes[0], _vertex_prototypes[1]),
        ),
        _LinearCageVertex.init_at_center(
            id=5,
            vertices=(_vertex_prototypes[0], _vertex_prototypes[2]),
        ),
        _LinearCageVertex.init_at_center(
            id=6,
            vertices=(_vertex_prototypes[0], _vertex_prototypes[3]),
        ),
        _LinearCageVertex.init_at_center(
            id=7,
            vertices=(_vertex_prototypes[1], _vertex_prototypes[2]),
        ),
        _LinearCageVertex.init_at_center(
            id=8,
            vertices=(_vertex_prototypes[1], _vertex_prototypes[3]),
        ),
        _LinearCageVertex.init_at_center(
            id=9,
            vertices=(_vertex_prototypes[2], _vertex_prototypes[3]),
        ),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[4]),
        Edge(1, _vertex_prototypes[0], _vertex_prototypes[5]),
        Edge(2, _vertex_prototypes[0], _vertex_prototypes[6]),
        Edge(3, _vertex_prototypes[1], _vertex_prototypes[4]),
        Edge(4, _vertex_prototypes[1], _vertex_prototypes[7]),
        Edge(5, _vertex_prototypes[1], _vertex_prototypes[8]),
        Edge(6, _vertex_prototypes[2], _vertex_prototypes[5]),
        Edge(7, _vertex_prototypes[2], _vertex_prototypes[7]),
        Edge(8, _vertex_prototypes[2], _vertex_prototypes[9]),
        Edge(9, _vertex_prototypes[3], _vertex_prototypes[6]),
        Edge(10, _vertex_prototypes[3], _vertex_prototypes[8]),
        Edge(11, _vertex_prototypes[3], _vertex_prototypes[9]),
    )

    _num_windows = 4
    _num_window_types = 1
