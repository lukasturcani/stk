from collections import Counter, defaultdict
import numpy as np
from functools import partial

from ..topology_graph import TopologyGraph, ConstructionState
from ...reactions import GenericReactionFactory


class _CageConstructionState(ConstructionState):
    def __init__(
        self,
        building_block_vertices,
        edges,
        scale,
        lattice_constants=None,
    ):
        super().__init__(
            building_block_vertices=building_block_vertices,
            edges=edges,
            scale=scale,
            lattice_constants=lattice_constants,
        )
        self._neighbor_positions = {}

    def clone(self):
        clone = super().clone()
        clone._neighbor_positions = {
            key: list(value)
            for key, value in self._neighbor_positions.items()
        }
        return clone

    def _with_neighbor_positions(self, neighbor_positions):
        for vertex_id, positions in neighbor_positions.items():
            self._neighbor_positions[vertex_id] = (
                self._neighbor_positions.get(vertex_id, [])
            )
            self._neighbor_positions[vertex_id].extend(
                np.array(position, dtype=np.float64)
                for position in positions
            )
        return self

    def with_neighbor_positions(self, neighbor_positions):
        return self.clone()._with_neighbor_positions(
            neighbor_positions=neighbor_positions,
        )

    def get_neighbor_positions(self, vertex_id):
        for position in self._neighbor_positions[vertex_id]:
            yield np.array(position, dtype=np.float64)


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

        """

        self._vertex_alignments = (
            dict(vertex_alignments)
            if vertex_alignments is not None
            else {}
        )

        super().__init__(
            vertices=tuple(
                vertex.with_aligner_edge(
                    aligner_edge=self._vertex_alignments.get(
                        vertex.get_id(),
                        0,
                    ),
                )
                for vertex in self._vertex_prototypes
            ),
            edges=self._edge_prototypes,
            reaction_factory=reaction_factory,
            construction_stages=tuple(
                partial(self._has_degree, degree)
                for degree
                in sorted(self._vertices_by_degree, reverse=True)
            ),
            num_processes=num_processes,
            edge_groups=None,
        )

    def _has_degree(self, degree, vertex):
        return vertex.get_id() in self._vertices_of_degree[degree]

    def get_building_block_vertices(self, building_blocks):
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
        for vertex in self._vertices:
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

    def _get_construction_state(self, building_block_vertices):
        return _CageConstructionState(
            building_block_vertices=building_block_vertices,
            edges=self._edges,
            scale=self._get_scale(building_block_vertices),
            lattice_constants=self._get_lattice_constants(),
        )

    def _after_placement_stage(
        self,
        state,
        vertices,
        edges,
        building_blocks,
        results,
    ):
        state = self._update_neighbor_positions(
            state=state,
            vertices=vertices,
            building_blocks=building_blocks,
            results=results,
        )
        return state.with_vertices(self._get_updated_vertices(state))

    def _update_neighbor_positions(
        self,
        state,
        vertices,
        edges,
        building_blocks,
        results,
    ):
        neighbor_positions = defaultdict(list)
        for vertex, vertex_edges, building_block, result in zip(
            vertices,
            edges,
            building_blocks,
            results,
        ):
            building_block = building_block.with_position_matrix(
                position_matrix=result.position_matrix,
            )
            edge_functional_groups = dict(zip(
                result.functional_group_edges.values(),
                result.functional_group_edges.keys(),
            ))
            for neighbor_id, edge_id in zip(
                self._get_neighbors(vertex, vertex_edges),
                (edge.get_id() for edge in vertex_edges),
            ):
                fg_id = edge_functional_groups[edge_id]
                functional_group = next(
                    building_block.get_functional_groups(fg_id)
                )
                neighbor_positions[neighbor_id].append(
                    building_block.get_centroid(
                        atom_ids=functional_group.get_placer_ids(),
                    )
                )

        return state.with_neighbor_positions(neighbor_positions)

    def _get_neighbors(self, vertex, vertex_edges):
        for edge in vertex_edges:
            yield (
                edge.get_vertex1_id()
                if vertex.get_id() != edge.get_vertex1_id()
                else edge.get_vertex2_id()
            )

    def _get_updated_vertices(self, state):
        for vertex_id in range(state.get_num_vertices()):
            neighbor_positions = tuple(state.get_neighbor_positions(
                vertex_id=vertex_id,
            ))
            if (
                len(neighbor_positions)
                == self._vertex_degrees[vertex_id]
            ):
                yield state.get_vertex(vertex_id).with_position(
                    position=(
                        sum(neighbor_positions)/len(neighbor_positions)
                    ),
                )
            else:
                yield state.get_vertex(vertex_id)

    def __repr__(self):
        vertex_alignments = (
            ', vertex_alignments={self._vertex_alignments}'
            if self._vertex_alignments
            else ''
        )
        return f'cage.{self.__class__.__name__}({vertex_alignments})'
