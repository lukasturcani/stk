"""
Edge Sorter
===========

"""

import numpy as np
from collections import abc

from ...edge import Edge
from .sorter import Sorter, IHasPosition


__all__ = (
    'EdgeSorter',
)


class EdgeSorter(Sorter):
    """
    Sorted edges according to their angle.

    """

    __slots__ = [
        '_items',
        '_reference',
        '_axis',
        '_edge_centroid',
    ]

    def __init__(
        self,
        edges: abc.Collection[Edge],
        aligner_edge: Edge,
        axis: np.ndarray,
    ) -> None:
        """
        Initialize an :class:`.EdgeSorter` instance.

        Parameters:
            edges:
                The edges to sort.

            aligner_edge:
                The edge in edges, used to calculate the reference
                vector.

            axis:
                Must be immutable. The axis used to determine the
                clockwise direction.

        """

        self._edge_centroid = edge_centroid = (
            sum(edge.get_position() for edge in edges) / len(edges)
        )
        super().__init__(
            items=edges,
            axis=axis,
            reference=aligner_edge.get_position() - edge_centroid,
        )

    def _get_vector(
        self,
        item: IHasPosition,
    ) -> np.ndarray:
        return item.get_position() - self._edge_centroid
