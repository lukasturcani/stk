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

from functools import partial
import numpy as np

from stk.utilities import flatten
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
        building_block_vertices,
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
        building_block_vertices : :class:`dict`
            Maps each :class:`.BuildingBlock` to be placed, to a
            :class:`tuple` of :class:`.Vertex` instances, on which
            it should be placed.

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

        scale = self._scale = self._get_scale(building_block_vertices)
        self._building_block_vertices = {
                building_block: tuple(
                    vertex.with_scale(scale) for vertex in vertices
                )
                for building_block, vertices
                in building_block_vertices.items()
        }
        self._edges = tuple(edge.with_scale(scale) for edge in edges)
        self._reaction_factory = reaction_factory
        if num_processes == 1:
            self._implementation = _Serial(
                stages=tuple(self._get_stages(construction_stages)),
            )
        else:
            self._implementation = _Parallel(
                stages=tuple(self._get_stages(construction_stages)),
                num_processes=num_processes,
            )

        if edge_groups is None:
            edge_groups = tuple(
                EdgeGroup((edge, )) for edge in self._edges
            )
        self._edge_groups = edge_groups

    def with_building_blocks(self, building_block_map):
        """
        Return a clone holding different building blocks.

        Parameters
        ----------
        building_block_map : :class:`dict`
            Maps a :class:`.BuildingBlock` in the current topology
            graph to the :class:`.BuildingBlock` which should replace
            it in the clone. If a building block should be not replaced
            in the clone, it can be omitted from the map.

        Returns
        -------
        :class:`.TopologyGraph`
            The clone.

        """

        ...

    def get_building_blocks(self):
        """
        Yield the building blocks.

        Building blocks are yielded in an order based on their
        position in the constructed molecule. For two equivalent
        topology graphs, but with different building blocks,
        equivalently positioned building blocks will be yielded at the
        same time.

        Yields
        ------
        :class:`.BuildingBlock`
            A building block of the topology graph.

        """

        yield from self._building_block_vertices

    def get_num_building_block(self, building_block):
        return len(
            self._building_block_vertices.get(building_block, [])
        )

    def _get_lattice_constants(self):
        return
        yield

    def construct(self):
        """
        Construct a :class:`.ConstructedMolecule`.

        Returns
        -------
        :class:`.ConstructionResult`
            The data describing the :class:`.ConstructedMolecule`.

        """

        state = self._get_construction_state()
        state = self._place_building_blocks(state)
        state = self._run_reactions(state)
        return ConstructionResult.init_from_construction_state(state)

    def _get_construction_state(self):
        return ConstructionState(
            building_block_vertices=self._building_block_vertices,
            edges=self._edges,
            lattice_constants=tuple(
                np.array(constant, dtype=np.float64)*self._scale
                for constant in self._get_lattice_constants()
            )
        )

    def _get_scale(self, building_block_vertices):
        raise NotImplementedError()

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
        vertices = flatten(self._building_block_vertices.values())
        for vertex in vertices:
            placed = False
            for i, stage in enumerate(construction_stages):
                if stage(vertex):
                    stages[i].append(vertex.get_id())
                    placed = True
                    break
            if not placed:
                stages[-1].append(vertex.get_id())
        yield from (stage for stage in stages if stage)

    def __str__(self):
        return repr(self)

    def __repr__(self):
        raise NotImplementedError()
