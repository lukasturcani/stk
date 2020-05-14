"""
Paddlewheel
===========

"""

from ..metal_complex import MetalComplex
from ..vertices import _MetalVertex, _BiDentateLigandVertex
from ...topology_graph import Edge


class Paddlewheel(MetalComplex):
    """
    Represents a metal complex topology graph.

    Metal building blocks with at least four functional groups are
    required for this topology.

    Ligand building blocks with two functional groups are required for
    this topology graph.

    When using a :class:`dict` for initialization, a
    :class:`.BuildingBlock` needs to be assigned to each of the
    following numbers:

        | metals: (0, 1)
        | ligands: (0, 1, 2, 3)

    See :class:`.MetalComplex` for more details and examples.

    """

    _metal_vertex_prototypes = (
        _MetalVertex(0, [0, 1, 0]),
        _MetalVertex(1, [0, -1, 0]),
    )
    _ligand_vertex_prototypes = (
        _BiDentateLigandVertex(2, [2, 0, 0]),
        _BiDentateLigandVertex(3, [0, 0, 2]),
        _BiDentateLigandVertex(4, [-2, 0, 0]),
        _BiDentateLigandVertex(5, [0, 0, -2]),
    )

    _edge_prototypes = (
        Edge(
            id=0,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[0],
            position=[0.1, 0.5, 0],
        ),
        Edge(
            id=1,
            vertex1=_metal_vertex_prototypes[1],
            vertex2=_ligand_vertex_prototypes[0],
            position=[0.1, -0.5, 0],
        ),

        Edge(
            id=2,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[1],
            position=[0, 0.5, 0.1],
        ),
        Edge(
            id=3,
            vertex1=_metal_vertex_prototypes[1],
            vertex2=_ligand_vertex_prototypes[1],
            position=[0, -0.5, 0.1],
        ),

        Edge(
            id=4,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[2],
            position=[-0.1, 0.5, 0],
        ),
        Edge(
            id=5,
            vertex1=_metal_vertex_prototypes[1],
            vertex2=_ligand_vertex_prototypes[2],
            position=[-0.1, -0.5, 0],
        ),

        Edge(
            id=6,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[3],
            position=[0, 0.5, -0.1],
        ),
        Edge(
            id=7,
            vertex1=_metal_vertex_prototypes[1],
            vertex2=_ligand_vertex_prototypes[3],
            position=[0, -0.5, -0.1],
        ),
    )
