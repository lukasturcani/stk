from .utilities import _Placement


class Serial:
    """
    Holds serial implementations of topology graph methods.

    """

    def __init__(self, stages, after_placement_stage):
        """
        Initialize a :class:`._Parallel`.

        Parameters
        ----------
        stages : :class:`tuple`
            A :class:`tuple` of the form ``([v1, v2, v3], [v4, v5])``,
            where each nested :class:`list` holds the
            :class:`.Vertex` objects used for placement in a particular
            stage.

        after_placement_stage : :class:`callable`
            The function to run on the :class:`._ConstructionState`
            after a placement stage has been completed. See
            :meth:`.TopologyGraph._after_placement_stage`.

        """

        self._stages = stages
        self._after_placement_stage = after_placement_stage

    def _place_building_blocks(self, state):
        vertex_assignments = state.get_vertex_assignments()
        for stage in self._stages:
            vertices = tuple(state.get_vertices(stage))
            building_blocks = tuple(
                map(vertex_assignments.get, stage)
            )
            edges = tuple(state.get_vertex_edges(stage))
            placements = map(
                _Placement,
                vertices,
                edges,
                building_blocks,
            )
            placement_results = map(
                lambda placement: placement.get_result(),
                placements,
            )
            state = state.with_placement_results(
                building_blocks=building_blocks,
                placement_results=placement_results,
            )
            state = self._after_placement_stage(
                state=state,
                vertices=vertices,
                edges=edges,
                building_blocks=building_blocks,
                placement_results=placement_results,
            )
        return state
