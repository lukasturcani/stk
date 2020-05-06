"""
Metal Complex
=============

"""

from itertools import chain

from ..topology_graph import TopologyGraph
from ...reactions import GenericReactionFactory


class MetalComplex(TopologyGraph):
    """
    Represents a metal complex topology graph.

    Examples
    --------
    *Construction*

    *Leaving Unsubstitued Sites*

    """

    def __init__(
        self,
        metals,
        ligands,
        reaction_factory=GenericReactionFactory(),
        num_processes=1,
    ):
        """
        Initialize a :class:`.MetalComplex`.

        Parameters
        ----------
        metals : :class:`dict`
            A :class:`dict` which maps the
            :class:`.BuildingBlock` instances to the ids of the
            vertices it should be placed on.

        ligands : :class:`dict`
            A :class:`dict` which maps the
            :class:`.BuildingBlock` instances to the ids of the
            vertices it should be placed on.

        reaction_factory : :class:`.ReactionFactory`, optional
            The reaction factory to use for creating bonds between
            building blocks.

        num_processes : :class:`int`, optional
            The number of parallel processes to create during
            :meth:`construct`.

        """

        metals = {
            metal: self._get_metal_vertices(ids)
            for metal, ids in metals.items()
        }
        ligands = {
            ligand: self._get_ligand_vertices(ids)
            for ligand, ids in ligands.items()
        }

        building_blocks = {**metals, **ligands}

        super().__init__(
            building_block_vertices={
                building_block: tuple(vertices)
                for building_block, vertices in building_blocks.items()
            },
            edges=self._edge_prototypes,
            reaction_factory=reaction_factory,
            construction_stages=(),
            num_processes=num_processes,
            edge_groups=None,
        )

    def _get_metal_vertices(self, vertex_ids):
        """
        Yield vertex prototypes.

        Parameters
        ----------
        vertex_ids : :class:`iterable` of :class:`int`
            The ids of the vertices to yield.

        Yields
        ------
        :class:`.Vertex`
            A vertex prototype of the topology graph.

        """

        if isinstance(vertex_ids, int):
            vertex_ids = (vertex_ids, )

        for vertex_id in vertex_ids:
            yield self._metal_vertex_prototypes[vertex_id]

    def _get_ligand_vertices(self, vertex_ids):
        """
        Yield vertex prototypes.

        Parameters
        ----------
        vertex_ids : :class:`iterable` of :class:`int`
            The ids of the vertices to yield.

        Yields
        ------
        :class:`.Vertex`
            A vertex prototype of the topology graph.

        """

        if isinstance(vertex_ids, int):
            vertex_ids = (vertex_ids, )

        for vertex_id in vertex_ids:
            yield self._ligand_vertex_prototypes[vertex_id]

    def clone(self):
        clone = super().clone()
        return clone

    def _get_scale(self, building_block_vertices):
        return 1

    def __repr__(self):
        return (
            f'metal_complex.{self.__class__.__name__}'
            f'()'
        )
