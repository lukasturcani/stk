"""
Serial Topology Graph
=====================

"""

from .utilities import _Placement


class _Serial:
    """
    Holds serial implementations of topology graph methods.

    """

    def __init__(self, stages):
        """
        Initialize a :class:`._Serial` instance.

        Parameters
        ----------
        stages : :class:`tuple`
            A :class:`tuple` of the form ``((v1, v2, v3), (v4, v5))``,
            where each nested :class:`tuple` holds the
            :class:`.Vertex` objects used for placement in a particular
            stage.

        """

        self._stages = stages

    def _place_building_blocks(self, state):
        for stage in self._stages:
            vertices = tuple(state.get_vertices(stage))
            building_blocks = tuple(
                map(state.get_building_block, stage)
            )
            edges = tuple(map(state.get_edges, stage))
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
                vertices=vertices,
                edges=edges,
                building_blocks=building_blocks,
                results=placement_results,
            )
        return state

    def get_num_stages(self):
        """
        Get the number of placement stages.

        Returns
        -------
        :class:`int`
            The number of placement stages.

        """

        return len(self._stages)
