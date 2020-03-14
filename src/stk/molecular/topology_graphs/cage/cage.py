from collections import Counter, defaultdict
from functools import partial
import numpy as np

from .cage_construction_state import _CageConstructionState
from ..topology_graph import TopologyGraph
from ...reactions import GenericReactionFactory


class Cage(TopologyGraph):
    """
    Represents a cage topology graph.

    Cage topologies are added by creating a subclass which defines the
    :attr:`_vertex_prototypes` and :attr:`_edge_prototypes` class
    attributes.

    Examples
    --------
    :class:`Cage` instances can be made without supplying
    additional arguments (using :class:`.FourPlusSix` as an example)

    .. code-block:: python

        import stk

        bb1 = stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()])
        bb2 = stk.BuildingBlock(
            smiles='O=CC(C=O)C=O',
            functional_groups=[stk.AldehydeFactory()],
        )
        cage1 = stk.ConstructedMolecule(
            building_blocks=(bb1, bb2),
            topology_graph=stk.cage.FourPlusSix()
        )

    Different structural isomers of cages can be made by using the
    `vertex_alignments` optional parameter

    .. code-block:: python

        tetrahedron = stk.cage.FourPlusSix(
            vertex_alignments={0: 1, 1: 1, 2: 2}
        )
        cage2 = stk.ConstructedMolecule(
            building_blocks=[bb1, bb2],
            topology_graph=tetrahedron
        )

    The parameter maps the id of a vertex to a number
    between 0 (inclusive) and the number of edges the vertex is
    connected to (exclusive). So a vertex connected to three edges
    can be mapped to ``0``, ``1`` or ``2``.

    By changing which edge each vertex is aligned with, a different
    structural isomer of the cage can be formed.

    You can also build cages with multiple building blocks, but you
    have to assign each building block to a vertex with
    `building_block_vertices`.

    bb1 = stk.BuildingBlock('O=CC(C=O)C=O', [stk.AldehydeFactory()])
    bb2 = stk.BuildingBlock(
        smiles='O=CC(Cl)(C=O)C=O',
        functional_groups=[stk.AldehydeFactory()],
    )
    bb3 = stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()])
    bb4 = stk.BuildingBlock('NCC(Cl)N', [stk.PrimaryAminoFactory()])
    bb5 = stk.BuildingBlock('NCCCCN', [stk.PrimaryAminoFactory()])

    tetrahedron = stk.cage.FourPlusSix()
    cage = stk.ConstructedMolecule(
        building_blocks=[bb1, bb2, bb3, bb4, bb5],
        topology_graph=tetrahedron,
        building_block_vertices={
            bb1: tetrahedron.get_vertices(range(2)),
            bb2: tetrahedron.get_vertices((2, 3)),
            bb3: tetrahedron.get_vertices(4),
            bb4: tetrahedron.get_vertices(5),
            bb5: tuple(tetrahedron.get_vertices())[6:],
        }
    )

    """

    def __init_subclass__(cls, **kwargs):
        cls._vertex_degrees = Counter(
            vertex_id
            for edge in cls._edge_prototypes
            for vertex_id in (
                edge.get_vertex1_id(),
                edge.get_vertex2_id(),
            )
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
        vertex_alignments : :class:`dict`, optional
            A mapping from the :attr:`.Vertex.id` of a :class:`.Vertex`
            :attr:`vertices` to an :class:`.Edge` connected to it.
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

        if isinstance(building_blocks, dict):
            building_blocks = {
                building_block: map(self._get_vertex, ids)
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
        with_aligner = partial(self._with_aligner, vertex_alignments)
        building_block_vertices = {
            building_block: tuple(map(with_aligner, vertices))
            for building_block, vertices in building_blocks.items()
        }
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

    def _get_vertex(self, vertex_id):
        return self._vertex_prototypes[vertex_id]

    @staticmethod
    def _with_aligner(vertex_alignments, vertex):
        """
        Return a clone of `vertex` with the aligner edge set.

        """

        return vertex.with_aligner_edge(
            aligner_edge=vertex_alignments.get(vertex.get_id(), 0),
        )

    def _has_degree(self, degree, vertex):
        return vertex.get_id() in self._vertices_of_degree[degree]

    def _get_building_block_vertices(self, building_blocks):
        bb_by_degree = {}
        for bb in building_blocks:
            if bb.get_num_functional_groups() in bb_by_degree:
                raise ValueError(
                    'If there are multiple building blocks with the '
                    'same number of functional groups, '
                    'building_block_vertices must be set explicitly.'
                )
            bb_by_degree[bb.get_num_functional_groups()] = bb

        building_block_vertices = {}
        for vertex in self._vertex_prototypes:
            bb = bb_by_degree[self._vertex_degrees[vertex.get_id()]]
            building_block_vertices[bb] = (
                building_block_vertices.get(bb, [])
            )
            building_block_vertices[bb].append(vertex)
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
            num_placement_stages=len(self._implementation._stages),
            vertex_degrees=self._vertex_degrees,
            lattice_constants=tuple(
                np.array(constant, dtype=np.float64)*self._scale
                for constant in self._get_lattice_constants()
            )
        )

    def __repr__(self):
        vertex_alignments = (
            ', vertex_alignments={self._vertex_alignments}'
            if self._vertex_alignments
            else ''
        )
        return f'cage.{self.__class__.__name__}({vertex_alignments})'
