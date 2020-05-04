"""
Metal Complex
=============

"""

from itertools import chain

from .vertices import (
    _MetalVertex,
    _MonoDentateLigandVertex,
    _BiDentateLigandVertex,
)
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
        building_blocks : :class:`iterable` or :class:`dict`
            Can be a :class:`iterable` of :class:`.BuildingBlock`
            instances, which should be placed on the topology graph.
            Can also be a :class:`dict` which maps the
            :class:`.BuildingBlock` instances to the ids of the
            vertices it should be placed on. A :class:`dict` is
            required when there are multiple building blocks with the
            same number of functional groups, because in this case
            the desired placement is ambiguous.

        vertex_alignments : :class:`dict`, optional
            A mapping from the id of a :class:`.Vertex`
            to an :class:`.Edge` connected to it.
            The :class:`.Edge` is used to align the first
            :class:`.FunctionalGroup` of a :class:`.BuildingBlock`
            placed on that vertex. Only vertices which need to have
            their default edge changed need to be present in the
            :class:`dict`. If ``None`` then the default edge is used
            for each vertex. Changing which :class:`.Edge` is used will
            mean that the topology graph represents different
            structural isomers. The edge is referred to by a number
            between ``0`` (inclusive) and the number of edges the
            vertex is connected to (exclusive).

        reaction_factory : :class:`.ReactionFactory`, optional
            The reaction factory to use for creating bonds between
            building blocks.

        num_processes : :class:`int`, optional
            The number of parallel processes to create during
            :meth:`construct`.

        Raises
        ------
        :class:`ValueError`
            If the there are multiple building blocks with the
            same number of functional_groups in `building_blocks`,
            and they are not explicitly assigned to vertices. The
            desired placement of building blocks is ambiguous in
            this case.

        """

        metals = {
            metal: self._get_metal_vertices(ids)
            for metal, ids in metals.items()
        }
        ligands = {
            ligand: self._get_ligand_vertices(ids)
            for ligand, ids in ligands.items()
        }
        print(metals, ligands)

        building_blocks_types = set(chain(
            metals.keys(), ligands.keys()
        ))
        print(building_blocks_types)
        building_blocks = {i: [] for i in building_blocks_types}
        for bb in building_blocks_types:
            if bb in metals:
                for v in metals[bb]:
                    building_blocks[bb].append(v)
            if bb in ligands:
                for v in ligands[bb]:
                    building_blocks[bb].append(v)

        print(building_blocks)

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

    # def _run_reactions(self, state):
    #     return state

    def _get_scale(self, building_block_vertices):
        return 3

    def __repr__(self):
        vertex_alignments = (
            f', vertex_alignments={self._vertex_alignments}'
            if self._vertex_alignments
            else ''
        )
        return (
            f'metal_complex.{self.__class__.__name__}'
            f'({vertex_alignments})'
        )
