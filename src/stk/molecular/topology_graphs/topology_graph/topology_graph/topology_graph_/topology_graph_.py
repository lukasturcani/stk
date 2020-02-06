import numpy as np
import pathos
from collections import namedtuple

from .implementations import _Parallel, _Serial
from ._construction_state import _ConstructionState
from ...utilities import vector_angle


class TopologyGraph_:
    """
    A partial implementation of :class:`.TopologyGraph`.

    """

    def __init__(
        self,
        vertex_data,
        edge_data,
        reaction_factory,
        construction_stages,
        num_processes,
    ):
        """
        Initialize an instance of :class:`.TopologyGraph`.

        Parameters
        ----------
        vertex_data : :class:`tuple` of :class:`.VertexData`
            The data for vertices which make up the graph.

        edge_data : :class:`tuple` of :class:`.EdgeData`
            The date for edges which make up the graph.

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

        """

        vertex_data = tuple(self._with_ids(vertex_data))
        edge_data = tuple(self._with_ids(edge_data))
        self._vertices = tuple(
            data.get_vertex() for data in vertex_data
        )
        self._edges = tuple(
            data.get_edge() for data in edge_data
        )
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

    def _get_construction_state(self, vertex_assignments):
        return _ConstructionState(


    def _after_placement_stage(
        self,
        state,
        vertices,
        edges,
        building_blocks,
        position_matrices,
        functional_group_to_edge_maps,
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

        position_matrices : :class:`tuple` of :class:`numpy.ndarray`
            For each building block in `building_blocks`, the position
            matrix after it has been placed by its vertex.

        functional_group_to_edge_maps : :class:`tuple` of :class:`dict`
            For each building block in `building_blocks`, the
            mapping of functional group to the edge id it is mapped to.

        Returns
        -------
        :class:`.ConstructionState`
            The new construction state.

        """

        return state

    @staticmethod
    def _with_ids(objects):
        for id, object in enumerate(objects):
            yield object.with_id(id)

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

    def _place_building_blocks(self, state):
        return self._implementation._place_building_blocks(state)

    def _run_reactions(self, state):
        reactions = map(
            self._reaction_factory.get_reaction,
            state.get_edge_functional_groups.values(),
        )
        results = map(
            lambda reaction: reaction.get_result(),
            reactions,
        )
        return state.with_reaction_results(results)
