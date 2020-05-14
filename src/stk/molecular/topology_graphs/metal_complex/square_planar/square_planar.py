"""
Square Planar
=============

"""

from ..metal_complex import MetalComplex
from ..vertices import _MetalVertex, _MonoDentateLigandVertex
from ...topology_graph import Edge


class SquarePlanar(MetalComplex):
    """
    Represents a square planar metal complex topology graph.

    Metal building blocks with at least four functional groups are
    required for this topology graph.

    Ligand building blocks with one functional group are required
    for this topology graph.

    When using a :class:`dict` for initialization, a
    :class:`.BuildingBlock` needs to be assigned to each of the
    following numbers:

        | metals: (0, )
        | ligands: (0, 1, 2, 3)

    See :class:`.MetalComplex` for more details and examples.

    """

    _metal_vertex_prototypes = (
        _MetalVertex(0, [0, 0, 0]),
    )
    _ligand_vertex_prototypes = (
        _MonoDentateLigandVertex(1, [2.5, 0, 0]),
        _MonoDentateLigandVertex(2, [0, 2.5, 0]),
        _MonoDentateLigandVertex(3, [-2.5, 0, 0]),
        _MonoDentateLigandVertex(4, [0, -2.5, 0]),
    )

    _edge_prototypes = (
        Edge(
            id=0,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[0],
        ),
        Edge(
            id=1,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[1],
        ),
        Edge(
            id=2,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[2],
        ),
        Edge(
            id=3,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[3],
        ),
    )
