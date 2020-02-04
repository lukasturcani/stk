import pathos


def _place_and_map(vertex, edges, building_block):
    """
    Place a `building_block` and map its functional groups to `edges`.

    Parameters
    ----------
    vertex : :class:`.Vertex`
        The vertex used to place the `building_block`.

    edges : :class:`tuple` of :class:`.Edge`
        The edges connected to the `vertex`.

    building_block : :class:`.BuildingBlock`
        The building block which is to be placed on the `vertex` and
        have its functional groups mapped to `edges`.

    Returns
    -------
    :class:`tuple`
        The position matrix of the `building_block` after it has
        been placed on the vertex is the first element of the
        :class:`tuple` and the mapping of functional groups to edge
        ids is the second element of the :class:`tuple`.

    """

    position_matrix = vertex.place_building_block(building_block)
    building_block = building_block.with_position_matrix(
        position_matrix=position_matrix,
    )
    functional_group_edges = vertex.map_functional_groups_to_edges(
        building_block=building_block,
        edges=edges,
    )
    return position_matrix, functional_group_edges


class _Parallel:
    """
    Holds parallel implementations of topology graph methods.

    """

    def __init__(self, stages, num_processes, after_placement_stage):
        """
        Initialize a :class:`._Parallel`.

        Parameters
        ----------
        stages : :class:`tuple`
            A :class:`tuple` of the form ``([v1, v2, v3], [v4, v5])``,
            where each nested :class:`list` holds the
            :class:`.Vertex` objects used for placement in a particular
            stage.

        num_processes : :class:`int`
            The number of parallel processes to spawn.

        after_placement_stage : :class:`callable`
            The function to run on the :class:`._ConstructionState`
            after a placement stage has been completed. See
            :meth:`.TopologyGraph._after_placement_stage`.

        """

        self._stages = stages
        self._num_processes = num_processes
        self._after_placement_stage = after_placement_stage

    def _place_building_blocks(self, state):
        vertex_assignments = state.get_vertex_assignments()

        with pathos.pools.ProcessPool(self._num_processes) as pool:
            for stage in self._stages:
                vertices = tuple(state.get_vertices(stage))
                building_blocks = tuple(
                    map(vertex_assignments.get, stage)
                )
                edges = tuple(state.get_vertex_edges(stage))
                position_matrices, maps = zip(*pool.map(
                    _place_and_map,
                    vertices,
                    edges,
                    building_blocks,
                ))
                state = state.with_building_blocks(
                    building_blocks=building_blocks,
                    position_matrices=position_matrices,
                    functional_group_to_edge_maps=maps,
                )
                state = self._after_placement_stage(
                    state=state,
                    vertices=vertices,
                    edges=edges,
                    building_blocks=building_blocks,
                    position_matrices=position_matrices,
                    maps=maps,
                )
        return state
