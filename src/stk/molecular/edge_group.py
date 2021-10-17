"""
Edge Group
==========

"""


from collections import abc
from .edge import Edge

__all__ = (
    'EdgeGroup',
)


class EdgeGroup:
    """
    Represents a group of edges.

    All edges in a group, and functional groups held by them, should be
    used together in a single :class:`.Reaction`.

    """

    def __init__(
        self,
        edges: abc.Iterable[Edge],
    ) -> None:
        """
        Initialize an :class:`.EdgeGroup` instance.

        Parameters:

            edges:
                The edges in the group.

        """

        self._edge_ids = tuple(edge.get_id() for edge in edges)

    def get_edge_ids(self) -> abc.Iterator[int]:
        """
        Yield the ids of edges in the group.

        Yields:

            The id of an edge.

        """

        yield from self._edge_ids

    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}(edge_ids={self._edge_ids})'
