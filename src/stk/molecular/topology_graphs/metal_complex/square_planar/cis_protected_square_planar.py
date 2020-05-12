"""
Cis Protected Square Planar
===========================

"""

from ..metal_complex import MetalComplex
from ..vertices import _MetalVertex, _BiDentateLigandVertex
from ...topology_graph import Edge


class CisProtectedSquarePlanar(MetalComplex):
    """
    Represents a square planar metal complex topology graph.

    Metal building blocks with at least two functional groups are
    required for this topology graph.

    Linker building blocks with two functional groups are required
    for this topology graph.

    .. code-block:: python

        # :class:`.BuildingBlock` placements:
        metals: (0, )
        ligands: (0, )

    See :class:`.MetalComplex` for more details and examples.

    """

    _metal_vertex_prototypes = (
        _MetalVertex(0, [0, 0, 0]),
    )
    _ligand_vertex_prototypes = (
        _BiDentateLigandVertex(1, [2.5, 2.5, 0]),
    )

    _edge_prototypes = (
        Edge(
            id=0,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[0],
            position=[2.5, 0, 0],
        ),
        Edge(
            id=1,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[0],
            position=[0, 2.5, 0],
        ),
    )
