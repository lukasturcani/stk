"""
Edge Sorter
===========

"""

from .sorter import _Sorter


class _EdgeSorter(_Sorter):
    """
    Sorted edges according to their angle.

    """

    __slots__ = [
        '_items',
        '_reference',
        '_axis',
        '_edge_centroid',
    ]

    def __init__(self, edges, aligner_edge, axis):
        """
        Initialize an :class:`._EdgeSorter` instance.

        Parameters
        ----------
        edges : :class:`iterable` of :class:`.Edge`
            The edges to sort.

        aligner_edge : :class:`.Edge`
            The edge in edges, used to calculate the reference vector.

        axis : :class:`numpy.ndarray`
            Must be immutable. The axis used to determine the clockwise
            direction.

        """

        self._edge_centroid = edge_centroid = (
            sum(edge.get_position() for edge in edges) / len(edges)
        )
        super().__init__(
            items=edges,
            axis=axis,
            reference=aligner_edge.get_position() - edge_centroid,
        )

    def _get_vector(self, item):
        return item.get_position() - self._edge_centroid
