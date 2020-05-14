"""
Porphyrin
=========

"""

from ..metal_complex import MetalComplex
from ..vertices import _MetalVertex
from ...topology_graph import Edge


class Porphyrin(MetalComplex):
    """
    Represents a metal complex topology graph.

    Metal building blocks with at least four functional groups are
    required for this topology.

    Ligand building blocks with two functional groups are required for
    this topology graph.

    When using a :class:`dict` for initialization, a
    :class:`.BuildingBlock` needs to be assigned to each of the
    following numbers:

        | metals: (0, )
        | ligands: (0, )

    See :class:`.MetalComplex` for more details and examples.

    """

    _metal_vertex_prototypes = (_MetalVertex(0, [0, 0, 0]), )
    _ligand_vertex_prototypes = (_MetalVertex(1, [0, 0, 0]), )

    _edge_prototypes = (
        Edge(
            id=0,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[0],
            position=[0.1, 0, 0],
        ),
        Edge(
            id=1,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[0],
            position=[0, 0.1, 0],
        ),
        Edge(
            id=2,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[0],
            position=[-0.1, 0, 0],
        ),
        Edge(
            id=3,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[0],
            position=[0, -0.1, 0],
        ),
    )
