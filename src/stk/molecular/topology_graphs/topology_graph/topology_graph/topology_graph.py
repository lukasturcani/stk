"""
Adding Topology Graphs
======================

To add a new topology graph a new subclass of :class:`.TopologyGraph`
must be added, which implements its virtual methods. Similarly, new
subclasses of :class:`.VertexData`` and :class:`.Vertex` must also be
made and their virtual methods implemented. When the new subclass of
:class:`.TopologyGraph` is initialized, it must create instances of the
:class:`.VertexData`
subclass, together with :class:`.EdgeData` instances. Once your
topology graph has the vertex and edge data it wants, simply run
:meth:`.TopologyGraph.__init__` and you're done.

When creating a :class:`.VertexData` subclass,
:meth:`~.VertexData.get_vertex` needs to be implemented such that it
returns an instance of your :class:`.Vertex` subclass. If you need to
define a new :meth:`__init__` method for either subclass, you will
also need to implement :meth:`clone` for it.

The :class:`.TopologyGraph` subclass can also create
`construction_stages` if parallel construction needs to be broken down
into separate stages. However,
if this is not the case, then an empty :class:`tuple` can simply be
passed.

Why is both :class:`.VertexData` and :class:`.Vertex` needed?
-------------------------------------------------------------

At first, it may appear that having both :class:`.VertexData` and
:class:`.Vertex` is an unnecessary inconvenience, as when you create
a new :class:`.TopologyGraph` subclass you have to subclass both of
these classes rather than just :class:`.Vertex`. The answer is
related to how these two classes reference other objects in the
:class:`.TopologyGraph`.

:class:`.VertexData` and :class:`.EdgeData` objects keep pointers
to each other in the :attr:`~.VertexData.edges` and
:attr:`~.EdgeData.vertices`. This is extremely convenient for
defining a :class:`.TopologyGraph` because its components can directly
reference each other. However, it poses a significant issue for
serialization. Topology graphs may be highly-cyclic structures
and are therefore they may not possible to serialize with off-the-shelf
serialization tools like :mod:`pickle` or :mod:`dill`. However,
serialization is necessary and fundamental for allowing
parallelization of :class:`.TopologyGraph` construction. The
vertices and edges of the graph have to be serialized and sent to
other cores so that they can place and connect building blocks in
parallel. As a  result, :class:`.VertexData` exists to allow a
convenient definition of a :class:`TopologyGraph`, while
:class:`.Vertex` exists to provide a serializable representation of it.
:class:`.Verex` and :class:`Edge` do not reference other objects
directly, instead they refer to them by their :attr:`.Vertex.id`,
which is used to get an index into :attr:`.TopologyGraph.vertices`
and :attr:`.TopologyGraph.edges`.

:meth:`.VertexData.get_vertex` is used to convert :class:`.VertexData`
into its :class:`.Vertex` counterpart.

"""

from collections import defaultdict
from functools import partial
import numpy as np

from ..construction_result import ConstructionResult
from ..construction_state import ConstructionState
from ..edge_group import EdgeGroup
from .implementations import _Parallel, _Serial


class TopologyGraph:
    """
    Represents topology graphs of :class:`.ConstructedMolecule`.

    The topology graph is an abstract representation of a constructed
    molecule. The vertices indicate where building blocks are placed
    and the edges indicate which building blocks have bonds formed
    between them by the construction process.

    Vertices are responsible for placing the building block molecules.
    By initializing the vertices with different parameters, you can
    alter how they position the building block molecules and therefore
    allow the user to easily specify a different structural isomer.

    Once a building block is placed on a vertex, the functional groups
    on the building block must be mapped to the different edges
    connected to the vertex. The number of functional groups in the
    building block must match the number of edges connected to the
    vertex.

    Once the functional groups are mapped to edges, each edge
    represents a reaction between the functional groups mapped to it.
    Note that more than two functional groups can map to the same edge,
    for cases where you are dealing with something really exotic.
    A :class:`.Reaction` between functional groups is selected based
    on the edges mapped to the edge. A :class:`.Reaction` will
    generally create bonds between the atoms of the functional groups.
    After this you will end up with a :class:`.ConstructedMolecule`.

    """

    def __init__(
        self,
        vertices,
        edges,
        reaction_factory,
        construction_stages,
        num_processes,
        edge_groups=None,
    ):
        """
        Initialize an instance of :class:`.TopologyGraph`.

        Parameters
        ----------
        vertices : :class:`tuple` of :class:`.VertexData`
            The vertices which make up the graph.

        edges : :class:`tuple` of :class:`.EdgeData`
            The edges which make up the graph.

        reaction_factory : :class:`.ReactionFactory`
            Used to pick which :class:`.Reaction` is used, given
            the functional groups on a topology graph edge.

        construction_stages : :class:`tuple` of :class:`callable`
            A collection of callables, each of which takes a
            :class:`.Vertex` and returns ``True`` or ``False``.
            If the first :class:`callable` is applied to a  vertex in
            `vertices`, that vertex is is part of the first
            construction stage. The second :class:`callable` is then
            applied to all vertices not in the first stage and those
            which return ``True`` belong to the second stage and
            so on.

            Vertices which belong to the same construction stage
            all place building blocks together in parallel, before
            placement is done by any vertices which are part of a later
            stage. This breaks down parallel construction into
            serial stages if synchronization between stages is needed.

            If the topology graph is performing construction serially,
            then all vertices which belong to an earlier stage will
            place their building block before those at a later stage.

        num_processes : :class:`int`
            The number of parallel processes to create during
            :meth:`construct`.

        edge_groups : :class:`tuple` of :class:`.EdgeGroup`, optional
            The edge groups of the topology graph, if ``None`` every
            :class:`.Edge` is in its own edge group.

        """

        self._vertices = vertices
        self._edges = edges
        self._vertex_edges = self._get_vertex_egdes(self._edges)
        self._reaction_factory = reaction_factory
        if num_processes == 1:
            self._implementation = _Serial(
                stages=self._get_stages(construction_stages),
                after_placement_stage=self._after_placement_stage,
            )
        else:
            self._implementation = _Parallel(
                stages=self._get_stages(construction_stages),
                num_processes=num_processes,
                after_placement_stage=self._after_placement_stage,
            )

        if edge_groups is None:
            edge_groups = tuple(
                EdgeGroup((edge, )) for edge in self._edges
            )
        self._edge_groups = edge_groups

    def _get_vertex_edges(self, edges):
        vertex_edges = defaultdict(list)
        for edge in edges:
            vertex_ids = (edge.get_vertex1_id(), edge.get_vertex2_id())

            if edge.is_periodic():
                for vertex_id in vertex_ids:
                    periodic_edge = self._get_periodic_edge(
                        edge=edge,
                        reference=vertex_id,
                    )
                    vertex_edges[vertex_id].append(periodic_edge)
            else:
                for vertex_id in vertex_ids:
                    vertex_edges[vertex_id].append(edge)

    def _get_periodic_edge(self, edge, reference):
        vertex1 = self._vertices[reference]
        vertex2 = self._vertices[
            edge.get_vertex1_id()
            if reference == edge.get_vertex2_id()
            else edge.get_vertex2_id()
        ]

        direction = 1 if reference == edge.get_vertex1_id() else -1
        periodicity = np.array(edge.get_periodicity())
        end_cell = (
            vertex1.get_cell() + direction*periodicity
        )

        cell_shift = end_cell - vertex2.get_cell()

        shift = 0
        lattice_constants = self._get_lattice_constants()
        for axis_shift, constant in zip(cell_shift, lattice_constants):
            shift += axis_shift*constant

        position = (
            (vertex2.get_position()+shift+vertex1.get_position()) / 2
        )
        return edge.with_position(position)

    def _get_lattice_constants(self):
        raise NotImplementedError()

    def construct(self, building_block_vertices):
        """
        Construct a :class:`.ConstructedMolecule`.

        Parameters
        ----------
        building_block_vertices : :class:`dict`
            Maps :class:`.BuildingBlock` instances to a
            :class:`tuple` holding :class:`.Vertex` instances.
            The mapping specifies which building block gets placed
            on which vertices of the graph.

        Returns
        -------
        :class:`.ConstructionResult`
            The data describing the :class:`.ConstructedMolecule`.

        """

        state = ConstructionState(
            building_block_vertices=building_block_vertices,
            edges=self._edges,
            vertex_edges=self._vertex_edges,
            scale=self._get_scale(building_block_vertices),
        )
        state = self._before_placement(state)
        state = self._place_building_blocks(state)
        state = self._before_reactions(state)
        state = self._run_reactions(state)
        state = self._clean_up(state)
        return ConstructionResult.init_from_construction_state(state)

    def get_vertices(self, vertex_ids=None):
        """
        Yield the vertices of the graph.

        Parameters
        ----------
        vertex_ids : :class:`iterable` of :class:`int`, optional
            The ids of vertices to yield. If ``None``, then all
            vertices are yielded.

        Yields
        ------
        :class:`.Vertex`
            A vertex of the graph.

        """

        if vertex_ids is None:
            vertex_ids = range(len(self._vertices))
        elif isinstance(vertex_ids, int):
            vertex_ids = (vertex_ids, )

        for id_ in vertex_ids:
            yield self._vertices[id_]

    def _get_scale(self, building_block_vertices):
        raise NotImplementedError()

    def _before_reactions(self, state):
        return state

    def _before_placement(self, state):
        return state

    def _place_building_blocks(self, state):
        """
        Place the building blocks onto the vertices.

        Parameters
        ----------
        state : :class:`._ConstructionState`
            Holds data necessary to construct the molecule.

        Returns
        -------
        :class:`._ConstructionState`
            The :class:`._ConstructionState` updated to account for
            the placed building blocks.

        """

        return self._implementation._place_building_blocks(state)

    def _run_reactions(self, state):
        get_reaction = partial(
            self._reaction_factory.get_reaction,
            state,
        )
        reactions = tuple(map(get_reaction, self._edge_groups))
        results = map(
            lambda reaction: reaction.get_result(),
            reactions,
        )
        return state.with_reaction_results(reactions, results)

    def _get_stages(self, construction_stages):
        stages = tuple(
            [] for i in range(len(construction_stages)+1)
        )
        for vertex in self._vertices:
            placed = False
            for i, stage in enumerate(construction_stages):
                if stage(vertex):
                    stages[i].append(vertex.get_id())
                    placed = True
                    break
            if not placed:
                stages[-1].append(vertex.get_id())
        return stages

    def _clean_up(self, state):
        return state

    def _after_placement_stage(
        self,
        state,
        vertices,
        edges,
        building_blocks,
        placement_results,
    ):
        """
        Perform `state` changes after a placement stage is done.

        This method should be overridden if the `state` needs to be
        modified after a placement stage of construction is complete.

        Parameters
        ----------
        state : :class:`._ConstructionState`
            The state of the construction process.

        vertices : :class:`tuple` of :class:`.Vertex`
            The vertices which were used in the last construction
            stage.

        edges : :class:`tuple`
            For each vertex in `vertices` a :class:`tuple` holding
            all connected edges. Has the form
            ``((e1, e2, e3), (e1, e4))``.

        building_blocks : :class:`tuple` of :class:`.BuildingBlock`
            The building blocks which were placed in the last
            construction stage. They are ordered to have the same
            index as the vertex in `vertices` which placed them.

        placement_results : :class:`.PlacementResult`
            The placement results of the previous placement stage.

        Returns
        -------
        :class:`.ConstructionState`
            The new construction state.

        """

        return state

    def __str__(self):
        return repr(self)

    def __repr__(self):
        raise NotImplementedError()
