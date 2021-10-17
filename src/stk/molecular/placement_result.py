"""
Placement Result
================

"""

from typing import NamedTuple
import numpy as np


__all__ = (
    'PlacementResult',
)


class PlacementResult(NamedTuple):
    """
    The result of a building block placement.

    Attributes:

        position_matrix:
            The position matrix of the building block after it has been
            placed on the vertex.

        functional_group_edges:
            Maps the id of a functional group to the id of the edge it
            is assigned to.

    """

    position_matrix: np.ndarray
    functional_group_edges: dict[int, int]
