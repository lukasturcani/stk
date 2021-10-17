"""
Serial Topology Graph
=====================

"""

from .utilities import Placement
from ...construction_state import ConstructionState

__all__ = (
    'Serial',
)


class Serial:
    """
    Holds serial implementations of topology graph methods.

    """

    def __init__(
        self,
        stages: tuple[tuple[int, ...], ...],
    ) -> None:
        """
        Initialize a :class:`.Serial` instance.

        Parameters:

            stages:
                A :class:`tuple` of the form
                ``((v1, v2, v3), (v4, v5))``, where each nested
                :class:`tuple` holds the id of the :class:`.Vertex`
                objects used for placement in a particular stage.

        """

        self._stages = stages

    def place_building_blocks(
        self,
        state: ConstructionState,
    ) -> ConstructionState:
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
            placement_results = map(
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
