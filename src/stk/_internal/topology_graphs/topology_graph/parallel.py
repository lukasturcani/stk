import pathos

from stk._internal.construction_state.construction_state import (
    ConstructionState,
)

from .utilities import _Placement


class _Parallel:
    """
    Holds parallel implementations of topology graph methods.
    """

    def __init__(
        self, stages: tuple[tuple[int, ...], ...], num_processes: int
    ) -> None:
        """
        Initialize a :class:`._Parallel` instance.

        Parameters:
            stages:
                A :class:`tuple` of the form ``((v1, v2, v3), (v4, v5))``,
                where each nested :class:`tuple` holds the
                :class:`.Vertex` objects used for placement in a particular
                stage.

            num_processes:
                The number of parallel processes to spawn.
        """
        self._stages = stages
        self._num_processes = num_processes

    def _place_building_blocks(
        self,
        state: ConstructionState,
    ) -> ConstructionState:
        with pathos.pools.ProcessPool(self._num_processes) as pool:
            for stage in self._stages:
                vertices = tuple(state.get_vertices(stage))
                building_blocks = tuple(map(state.get_building_block, stage))
                edges = tuple(map(state.get_edges, stage))
                placements = map(
                    _Placement,
                    vertices,
                    edges,
                    building_blocks,
                )
                placement_results = pool.map(
                    lambda placement: placement.get_result(),
                    placements,
                )
                state = state.with_placement_results(
                    vertices=vertices,
                    edges=edges,
                    building_blocks=building_blocks,
                    results=placement_results,
                )
        return state

    def get_num_stages(self) -> int:
        """
        Get the number of placement stages.

        Returns:
            The number of placement stages.
        """
        return len(self._stages)
