"""
Placement Result
================

"""

import numpy as np


__all__ = (
    'PlacementResult',
)


class PlacementResult:
    """
    The result of a building block placement.

    """

    __slots__ = (
        'position_matrix',
        'functional_group_edges',
    )

    def __init__(
        self,
        position_matrix: np.ndarray,
        functional_group_edges: dict[int, int],
    ) -> None:
        """
        Initialize a :class:`.PlacementResult`.

        Parameters:

            position_matrix:
                The position matrix of the building block after it has
                been placed on a vertex.

            functional_group_edges:
                Maps the id of a functional group to the id of the edge
                it is assigned to.
        """

        self.position_matrix = position_matrix
        self.functional_group_edges = functional_group_edges
