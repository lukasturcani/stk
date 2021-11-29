"""
Topology Graph
==============

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
    An abstract base class for topology graphs.

    It is responsible for the construction of molecules. To create
    a new topology graph, you want to subclass and implement this
    abstract base class.

    Notes
    -----
    *Adding New Topology Graphs*

    You might notice that some of the methods of this abstract base
    class are implemented. This is purely for convenience when
    implementing subclasses. The implemented public methods are
    simply default implementations, which can safely be ignored or
    overridden, when implementing subclasses. Any private methods are
    implementation details of these default implementations.

    Many classes, such as :class:`.Vertex`, :class:`.Edge`,
    :class:`.EdgeGroup` and :class:`.ConstructionState`, exist as
    implementation details of this default :class:`.TopologyGraph`
    implementation. You could ignore all of them, and define a new
    :meth:`.construct` method from scratch. In fact, your topology
    graph does not have to be represented as a graph at all. However,
    using the default implementation of :class:`.TopologyGraph` makes
    it significantly easier to implement a construction process. When
    using the default implementation of :class:`.TopologyGraph`, you
    mostly just need to implement a :class:`.Vertex` subclass, which
    is much easier than figuring out the whole construction process
    from scratch. In addition, you get benefits like parallel
    construction for free, as it is included in the default
    implementation.

    Typically, adding a new topology graph will involve implementing
    any pure virtual methods of :class:`.TopologyGraph`, in a new
    subclass, as well as implementing any pure virtual methods of
    :class:`.Vertex`, again in a new subclass. Combined, this is just a
    handful of simple methods to implement. Sometimes, rarely, you
    might also want to subclass :class:`.ConstructionState`, when you
    want to add additional hooks during construction, by extending
    the methods of this class. If you do this, make sure
    to override :meth:`._get_construction_state` to return your
    subclass of :class:`.ConstructionState`, rather than the base
    class, as is done by default. You can subclass and extend the
    methods of any class as you wish, but it would be unusual if this
    doesn't cover all your requirements.

    *The Default Implementation*

    The default implementation of :class:`.TopologyGraph` represents
    the constructed molecule through a graph. The vertices indicate
    where building blocks are placed and the edges indicate which
    building blocks have bonds formed between them by the construction
    process.

    :class:`.Vertex` instances are responsible for placing the building
    block molecules. By initializing the vertices with different
    parameters, you can alter how they position the building block
    molecules, and therefore allow the user to easily specify a
    different structural isomer.

    Once a building block is placed on a vertex, the functional groups
    on the building block must be mapped to the different edges
    connected to the vertex. The number of functional groups in the
    building block must match the number of edges connected to the
    vertex.

    Once the functional groups are mapped to edges, the edges are
    used to perform reactions on the building blocks. Edges are
    grouped in an :class:`.EdgeGroup`, and all functional groups
    present in the edge group are reacted together. Normally, unless
    you are doing something very exotic, an :class:`.EdgeGroup` will
    hold just one :class:`.Edge`, and the two functional groups on
    that edge will be reacted together through a single
    :class:`.Reaction`. This reaction will normally add the bonds which
    are required to form the joined-up constructed molecule, but note
    that it does not have to add any bonds at all. In addition, a
    :class:`.Reaction` can add and remove atoms from the constructed
    molecule. Which reaction is selected to join the functional groups
    depends on the :class:`.ReactionFactory` given to the
    :class:`.TopologyGraph` during initialization.

    Once this is done, you have a :class:`.ConstructedMolecule`.

    Examples
    --------
    *Subclass Implementation*

    The source code of subclasses, listed in
    :mod:`~.topology_graph.topology_graph.topology_graph`, can serve
    as good examples.

    *Changing the Building Blocks of a Topology Graph*

    To change the building blocks used by a topology graph you
    can use :meth:`.with_building_blocks` to get a clone of the
    topology graph holding the new building blocks

    .. testcode:: changing-the-building-blocks-of-a-topology-graph

        import stk

        bb1 = stk.BuildingBlock('BrCCBr', [stk.BromoFactory()])
        bb2 = stk.BuildingBlock('BrCCCBr', [stk.BromoFactory()])

        linear = stk.polymer.Linear(
            building_blocks=(bb1, bb2),
            repeating_unit='A',
            num_repeating_units=15,
        )

        bb3 = stk.BuildingBlock('BrCNCBr', [stk.BromoFactory()])
        # All bb1 instances are replaced by bb3, but bb2 remains
        # in place.
        clone = linear.with_building_blocks({
            bb1: bb3,
        })

    """

    def __init__(
        self,
        building_block_vertices,
        edges,
        reaction_factory,
        construction_stages,
        num_processes,
        optimizer,
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

        edges : :class:`tuple` of :class:`.Edge`
            The edges which make up the topology graph.

        reaction_factory : :class:`.ReactionFactory`
            Used to pick which :class:`.Reaction` is used on each
            :class:`.EdgeGroup` of the topology graph.

        construction_stages : :class:`tuple` of :class:`callable`
            A collection of callables, each of which takes a
            :class:`.Vertex` and returns ``True`` or ``False``.
            If the first :class:`callable` is applied to a  vertex in
            the topology graph, and the result is ``True``, that vertex
            is a part of the first construction stage. The second
            :class:`callable` is then applied to all vertices not in
            the first stage and those which return ``True`` belong to
            the second stage and so on.

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

        optimizer : :class:`.Optimizer`
            Used to optimize the structure of the constructed
            molecule.

        edge_groups : :class:`tuple` of :class:`.EdgeGroup`, optional
            The edge groups of the topology graph, if ``None``, every
            :class:`.Edge` is in its own edge group.

        """

        self._scale = scale = self._get_scale(building_block_vertices)

        def apply_scale(item):
            return item.with_scale(scale)

        self._building_block_vertices = {
            building_block: tuple(map(apply_scale, vertices))
            for building_block, vertices
            in building_block_vertices.items()
        }
        self._edges = tuple(map(apply_scale, edges))
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

        self._optimizer = optimizer

    def _with_building_blocks(self, building_block_map):
        """
        Modify the topology graph.

        """

        # The original scaling first needs to be removed, so that when
        # the scale is recalculated with the new building blocks, it
        # has the same starting geometry.
        def undo_scale(vertex):
            return vertex.with_scale(1/self._scale)

        building_block_vertices = {
            building_block_map.get(building_block, building_block):
                tuple(map(undo_scale, vertices))
            for building_block, vertices
            in self._building_block_vertices.items()
        }
        scale = self._get_scale(building_block_vertices)

        def scale_vertex(vertex):
            return vertex.with_scale(scale)

        self._building_block_vertices = {
            building_block: tuple(map(scale_vertex, vertices))
            for building_block, vertices
            in building_block_vertices.items()
        }

        def scale_edge(edge):
            # Remove the old scale and apply the new one.
            return edge.with_scale(scale/self._scale)

        self._edges = edges = tuple(map(scale_edge, self._edges))

        def get_new_edge(edge_id):
            return edges[edge_id]

        self._edge_groups = tuple(
            EdgeGroup(map(get_new_edge, edge_group.get_edge_ids()))
            for edge_group in self._edge_groups
        )
        self._scale = scale
        return self

    def with_building_blocks(self, building_block_map):
        """
        Return a clone holding different building blocks.

        Parameters
        ----------
        building_block_map : :class:`dict`
            Maps a building block in the current topology
            graph to the building block which should replace
            it in the clone. If a building block should be not replaced
            in the clone, it can be omitted from the map.

        Returns
        -------
        :class:`.TopologyGraph`
            The clone. Has the same type as the original topology
            graph.

        """

        return self.clone()._with_building_blocks(building_block_map)

    def clone(self):
        """
        Return a clone.

        Returns
        -------
        :class:`.TopologyGraph`
            The clone. Has the same type as the original topology
            graph.

        """

        clone = self.__class__.__new__(self.__class__)
        clone._scale = self._scale
        clone._building_block_vertices = dict(
            self._building_block_vertices
        )
        clone._edges = self._edges
        clone._reaction_factory = self._reaction_factory
        clone._implementation = self._implementation
        clone._optimizer = self._optimizer
        clone._edge_groups = self._edge_groups
        return clone

    def get_building_blocks(self):
        """
        Yield the building blocks.

        Building blocks are yielded in an order based on their
        position in the topology graph. For two equivalent
        topology graphs, but with different building blocks,
        equivalently positioned building blocks will be yielded at the
        same time.

        Yields
        ------
        :class:`.BuildingBlock`
            A building block of the topology graph.

        """

        vertex_building_blocks = {}
        num_vertices = 0
        for building_block, vertices in (
            self._building_block_vertices.items()
        ):
            for vertex in vertices:
                num_vertices += 1
                vertex_building_blocks[vertex.get_id()] = (
                    building_block
                )

        yielded = set()
        for vertex_id in range(num_vertices):
            building_block = vertex_building_blocks[vertex_id]
            if building_block not in yielded:
                yielded.add(building_block)
                yield building_block

    def get_num_building_block(self, building_block):
        """
        Get the number of times `building_block` is present.

        Parameters
        ----------
        building_block : :class:`.BuildingBlock`
            The building block whose frequency in the topology graph
            is desired.

        Returns
        -------
        :class:`int`
            The number of times `building_block` is present in the
            topology graph.

        """

        return len(
            self._building_block_vertices.get(building_block, [])
        )

    def _get_lattice_constants(self):
        """
        Yield the lattice constants of the topology graph.

        The a, b and c lattice constants are yielded, in that order.

        By default, this is an empty generator.

        Yields
        ------
        :class:`numpy.ndarray`
            A lattice constant.

        """

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
        state = self._optimizer.optimize(state)
        return self._get_construction_result(state)

    def _get_construction_result(self, state):
        """
        Get the result of the construction.

        Parameters
        ----------
        state : :class:`.ConstructionState`
            The state of the molecule being constructed.

        Returns
        -------
        :class:`.ConstructionResult`
            The data describing the :class:`.ConstructedMolecule`.

        """

        return ConstructionResult(state)

    def _get_construction_state(self):
        return ConstructionState(
            building_block_vertices=self._building_block_vertices,
            edges=self._edges,
            lattice_constants=tuple(
                np.array(constant, dtype=np.float64)*self._scale
                for constant in self._get_lattice_constants()
            ),
        )

    def _get_scale(self, building_block_vertices):
        """
        Get the scale, which should be applied to topology graph.

        The scale should be applied to the position of every vertex
        and edge of topology graph. This allows to graph to adjust
        based on the size of the building blocks.

        Parameters
        ----------
        building_block_vertices : :class:`dict`
            Maps every :class:`.BuildingBlock` of the topology graph,
            to a :class:`tuple` of the :class:`.Vertex` instances it
            is meant to be placed on.

        Returns
        -------
        :class:`float`
            The scale.

        """

        raise NotImplementedError()

    def _place_building_blocks(self, state):
        """
        Place the building blocks onto the vertices.

        Parameters
        ----------
        state : :class:`.ConstructionState`
            Holds data necessary to construct the molecule.

        Returns
        -------
        :class:`.ConstructionState`
            The new construction state, updated to account for the
            placed building blocks.

        """

        return self._implementation._place_building_blocks(state)

    def _run_reactions(self, state):
        """
        Perform the reactions on the building blocks.

        Parameters
        ----------
        state : :class:`.ConstructionState`
            The current state of the construction process.

        Returns
        -------
        :class:`.ConstructionState`
            The new construction state, updated to account for the
            reactions between building blocks.

        """

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
        """
        Yield the parallelizable stages of construction.

        Yields
        ------
        :class:`tuple` of :class:`.Vertex`
            Vertices, which can be placed in parallel.

        """

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
        yield from (tuple(stage) for stage in stages if stage)

    def __str__(self):
        return repr(self)

    def __repr__(self):
        raise NotImplementedError()
