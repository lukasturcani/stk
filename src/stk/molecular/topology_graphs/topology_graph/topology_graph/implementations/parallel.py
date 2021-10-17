"""
Parallel Topology Graph
=======================

"""

import pathos
import typing

from .utilities import Placement
from ...construction_state import ConstructionState
from ...vertex import Vertex

__all__ = (
    'Parallel',
)

_V = typing.TypeVar('_V', bound=Vertex)


class Parallel:
    """
    Holds parallel implementations of topology graph methods.

    """

    def __init__(
        self,
        stages: tuple[tuple[int, ...], ...],
        num_processes: int,
    ) -> None:
        """
        Initialize a :class:`.Parallel` instance.

        Parameters:

            stages:
                A :class:`tuple` of the form
                ``((v1, v2, v3), (v4, v5))``, where each nested
                :class:`tuple` holds the id of the :class:`.Vertex`
                objects used for placement in a particular stage.

            num_processes:
                The number of parallel processes to spawn.

        """

        self._stages = stages
        self._num_processes = num_processes

    def place_building_blocks(
        self,
        state: ConstructionState[Vertex],
    ) -> ConstructionState:
        with pathos.pools.ProcessPool(self._num_processes) as pool:
            for stage in self._stages:
                vertices = tuple(state.get_vertices(stage))
                building_blocks = tuple(
                    map(state.get_building_block, stage)
                )
                edges = tuple(map(state.get_edges, stage))
                placements = map(
                    Placement,
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
