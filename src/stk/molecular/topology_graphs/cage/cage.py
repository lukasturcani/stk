"""
Cage
====

.. toctree::
    :maxdepth: 2

    One Plus One <\
stk.molecular.topology_graphs.cage.three_plus_three.one_plus_one\
>
    Two Plus Two <\
stk.molecular.topology_graphs.cage.three_plus_three.two_plus_two\
>
    Two Plus Three <\
stk.molecular.topology_graphs.cage.two_plus_three.two_plus_three\
>
    Two Plus Four <\
stk.molecular.topology_graphs.cage.two_plus_four.two_plus_four\
>
    Three Plus Six <\
stk.molecular.topology_graphs.cage.two_plus_four.three_plus_six\
>
    Four Plus Four <\
stk.molecular.topology_graphs.cage.three_plus_three.four_plus_four\
>
    Four Plus Six <\
stk.molecular.topology_graphs.cage.two_plus_three.four_plus_six\
>
    Four Plus Six 2 <\
stk.molecular.topology_graphs.cage.two_plus_three.four_plus_six_2\
>
    Four Plus Eight <\
stk.molecular.topology_graphs.cage.two_plus_four.four_plus_eight\
>
    Five Plus Ten <\
stk.molecular.topology_graphs.cage.two_plus_four.five_plus_ten\
>
    Six Plus Eight <\
stk.molecular.topology_graphs.cage.three_plus_four.six_plus_eight\
>
    Six Plus Nine <\
stk.molecular.topology_graphs.cage.two_plus_three.six_plus_nine\
>
    Six Plus Twelve <\
stk.molecular.topology_graphs.cage.two_plus_four.six_plus_twelve\
>
    Eight Plus Twelve <\
stk.molecular.topology_graphs.cage.two_plus_three.eight_plus_twelve\
>
    Eight Plus Sixteen <\
stk.molecular.topology_graphs.cage.two_plus_four.eight_plus_sixteen\
>
    Ten Plus Twenty <\
stk.molecular.topology_graphs.cage.two_plus_four.ten_plus_twenty\
>
    Twelve Plus Thirty <\
stk.molecular.topology_graphs.cage.two_plus_five.twelve_plus_thirty\
>
    Twenty Plus Thirty <\
stk.molecular.topology_graphs.cage.two_plus_three.twenty_plus_thirty\
>

"""

from collections import Counter, defaultdict
from functools import partial

from .cage_construction_state import _CageConstructionState
from ..topology_graph import TopologyGraph
from ...reactions import GenericReactionFactory


class UnoccupiedVertexError(Exception):
    """
    When a cage vertex is not occupied by a building block.

    """

    pass


class OverlyOccupiedVertexError(Exception):
    """
    When a cage vertex is occupied by more than one building block.

    """

    pass


class Cage(TopologyGraph):
    """
    Represents a cage topology graph.

    Notes
    -----
    Cage topologies are added by creating a subclass, which defines the
    :attr:`_vertex_prototypes` and :attr:`_edge_prototypes` class
    attributes.

    Examples
    --------
    *Subclass Implementation*

    The source code of the subclasses, listed in :mod:`~.cage.cage`,
    can serve as good examples.

    *Basic Construction*

    :class:`.Cage` instances can be made by providing the building
    block molecules only (using :class:`.FourPlusSix` as an example)

    .. code-block:: python

        import stk

        bb1 = stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()])
        bb2 = stk.BuildingBlock(
            smiles='O=CC(C=O)C=O',
            functional_groups=[stk.AldehydeFactory()],
        )
        cage1 = stk.ConstructedMolecule(
            topology_graph=stk.cage.FourPlusSix((bb1, bb2)),
        )

    *Structural Isomer Construction*

    Different structural isomers of cages can be made by using the
    `vertex_alignments` optional parameter

    .. code-block:: python

        cage2 = stk.ConstructedMolecule(
            topology_graph=stk.cage.FourPlusSix(
                building_blocks=(bb1, bb2),
                vertex_alignments={0: 1, 1: 1, 2: 2},
            ),
        )

    The parameter maps the id of a vertex to a number
    between 0 (inclusive) and the number of edges the vertex is
    connected to (exclusive). So a vertex connected to three edges
    can be mapped to ``0``, ``1`` or ``2``.

    By changing which edge each vertex is aligned with, a different
    structural isomer of the cage can be formed.

    *Multi-Building Block Cage Construction*

    You can also build cages with multiple building blocks, but,
    if you have multiple building blocks with the same number
    of functional groups, you have to assign each building block to the
    vertex you want to place it on

    .. code-block:: python

        bb1 = stk.BuildingBlock(
            smiles='O=CC(C=O)C=O',
            functional_groups=[stk.AldehydeFactory()],
        )
        bb2 = stk.BuildingBlock(
            smiles='O=CC(Cl)(C=O)C=O',
            functional_groups=[stk.AldehydeFactory()],
        )
        bb3 = stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()])
        bb4 = stk.BuildingBlock(
            smiles='NCC(Cl)N',
            functional_groups=[stk.PrimaryAminoFactory()],
        )
        bb5 = stk.BuildingBlock('NCCCCN', [stk.PrimaryAminoFactory()])

        cage = stk.ConstructedMolecule(
            topology_graph=stk.cage.FourPlusSix(
                # building_blocks is now a dict, which maps building
                # blocks to the id of the vertices it should be placed
                # on. You can use ranges to specify the ids.
                building_blocks={
                    bb1: range(2),
                    bb2: (2, 3),
                    bb3: 4,
                    bb4: 5,
                    bb5: range(6, 10),
                },
            ),
        )

    You can combine this with the `vertex_alignments` parameter

    .. code-block:: python

        cage = stk.ConstructedMolecule(
            topology_graph=stk.cage.FourPlusSix(
                building_blocks={
                    bb1: range(2),
                    bb2: (2, 3),
                    bb3: 4,
                    bb4: 5,
                    bb5: range(6, 10),
                },
            ),
            vertex_alignments={0: 1, 1: 1, 2: 2},
        )

    """

    def __init_subclass__(cls, **kwargs):
        cls._vertex_degrees = Counter(
            vertex_id
            for edge in cls._edge_prototypes
            for vertex_id in edge.get_vertex_ids()
        )
        cls._vertices_of_degree = defaultdict(set)
        for vertex_id, degree in cls._vertex_degrees.items():
            cls._vertices_of_degree[degree].add(vertex_id)

    def __init__(
        self,
        building_blocks,
        vertex_alignments=None,
        reaction_factory=GenericReactionFactory(),
        num_processes=1,
    ):
        """
        Initialize a :class:`.Cage`.

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
        :class:`AssertionError`
            If the any building block does not have a
            valid number of functional groups.

        :class:`ValueError`
            If the there are multiple building blocks with the
            same number of functional_groups in `building_blocks`,
            and they are not explicitly assigned to vertices. The
            desired placement of building blocks is ambiguous in
            this case.

        :class:`~.cage.UnoccupiedVertexError`
            If a vertex of the cage topology graph does not have a
            building block placed on it.

        :class:`~.cage.OverlyOccupiedVertexError`
            If a vertex of the cage topology graph has more than one
            building block placed on it.

        """

        # Use tuple here because it prints nicely.
        allowed_degrees = tuple(self._vertices_of_degree.keys())
        if isinstance(building_blocks, dict):
            for building_block in building_blocks:
                assert (
                    building_block.get_num_functional_groups()
                    in self._vertices_of_degree.keys()
                ), (
                    'The number of functional groups in '
                    f'{building_block} needs to be one of '
                    f'{allowed_degrees}, but is '
                    'currently '
                    f'{building_block.get_num_functional_groups()}.'
                )
            building_blocks = {
                building_block: self._get_vertices(ids)
                for building_block, ids in building_blocks.items()
            }

        else:
            building_blocks = self._get_building_block_vertices(
                building_blocks=building_blocks,
            )

        self._vertex_alignments = vertex_alignments = (
            dict(vertex_alignments)
            if vertex_alignments is not None
            else {}
        )

        def with_aligner(vertex):
            return vertex.with_aligner_edge(
                aligner_edge=vertex_alignments.get(vertex.get_id(), 0),
            )

        building_block_vertices = {
            building_block: tuple(map(with_aligner, vertices))
            for building_block, vertices in building_blocks.items()
        }
        self._check_building_block_vertices(building_block_vertices)
        super().__init__(
            building_block_vertices=building_block_vertices,
            edges=self._edge_prototypes,
            reaction_factory=reaction_factory,
            construction_stages=tuple(
                partial(self._has_degree, degree)
                for degree
                in sorted(self._vertices_of_degree, reverse=True)
            ),
            num_processes=num_processes,
            edge_groups=None,
        )

    @classmethod
    def _check_building_block_vertices(cls, building_block_vertices):
        unassigned_ids = set(
            vertex.get_id() for vertex in cls._vertex_prototypes
        )
        assigned_ids = set()
        vertices = (
            vertex
            for vertices_ in building_block_vertices.values()
            for vertex in vertices_
        )
        for vertex in vertices:
            if vertex.get_id() in assigned_ids:
                raise OverlyOccupiedVertexError(
                    f'Vertex {vertex.get_id()} has multiple building '
                    'blocks placed on it.'
                )
            assigned_ids.add(vertex.get_id())
            unassigned_ids.remove(vertex.get_id())

        if unassigned_ids:
            raise UnoccupiedVertexError(
                'The following vertices are unoccupied '
                f'{unassigned_ids}.'
            )

    def clone(self):
        clone = super().clone()
        clone._vertex_alignments = dict(self._vertex_alignments)
        return clone

    def _get_vertices(self, vertex_ids):
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
            yield self._vertex_prototypes[vertex_id]

    def _has_degree(self, degree, vertex):
        """
        Check if `vertex` has a degree of `degree`.

        Parameters
        ----------
        degree : :class:`int`
            The degree in question.

        vertex : :class:`.Vertex`
            The vertex in question.

        Returns
        -------
        :class:`bool`
            ``True`` if `vertex` has a degree of `degree`.

        """

        return vertex.get_id() in self._vertices_of_degree[degree]

    def _get_building_block_vertices(self, building_blocks):
        """
        Map building blocks to the vertices of the graph.

        Parameters
        ----------
        building_blocks : :class:`iterable` of :class:`.BuildingBlock`
            The building blocks which need to be mapped to vertices.

        Returns
        -------
        :class:`dict`
            Maps each building block in `building_blocks` to a
            :class:`list` of :class:`.Vertex` instances it should be
            placed on.

        Raises
        ------
        :class:`AssertionError`
            If the any building block does not have a
            valid number of functional groups.

        :class:`ValueError`
            If there are multiple building blocks with the same number
            of functional groups.

        """

        # Use tuple here because it prints nicely.
        allowed_degrees = tuple(self._vertices_of_degree.keys())

        building_blocks_by_degree = {}
        for building_block in building_blocks:
            num_fgs = building_block.get_num_functional_groups()
            assert (
                num_fgs in self._vertices_of_degree.keys()
            ), (
                'The number of functional groups in '
                f'{building_block} needs to be one of '
                f'{allowed_degrees}, but is '
                'currently '
                f'{building_block.get_num_functional_groups()}.'
            )

            if num_fgs in building_blocks_by_degree:
                raise ValueError(
                    'If there are multiple building blocks with the '
                    'same number of functional groups, '
                    'building_block_vertices must be set explicitly.'
                )
            building_blocks_by_degree[num_fgs] = building_block

        building_block_vertices = {}
        for vertex in self._vertex_prototypes:
            vertex_degree = self._vertex_degrees[vertex.get_id()]
            building_block = building_blocks_by_degree[vertex_degree]
            building_block_vertices[building_block] = (
                building_block_vertices.get(building_block, [])
            )
            building_block_vertices[building_block].append(vertex)
        return building_block_vertices

    def _get_scale(self, building_block_vertices):
        return max(
            bb.get_maximum_diameter()
            for bb in building_block_vertices
        )

    def _get_construction_state(self):
        return _CageConstructionState(
            building_block_vertices=self._building_block_vertices,
            edges=self._edges,
            num_placement_stages=self._implementation.get_num_stages(),
            vertex_degrees=self._vertex_degrees,
        )

    def __repr__(self):
        vertex_alignments = (
            ', vertex_alignments={self._vertex_alignments}'
            if self._vertex_alignments
            else ''
        )
        return f'cage.{self.__class__.__name__}({vertex_alignments})'
