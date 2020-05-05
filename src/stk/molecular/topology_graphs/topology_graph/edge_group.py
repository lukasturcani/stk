"""
Edge Group
==========

"""


class EdgeGroup:
    """
    Represents a group of edges.

    All edges in a group, and functional groups held by them, should be
    used together in a single :class:`.Reaction`.

    """

    def __init__(self, edges):
        """
        Initialize an :class:`.EdgeGroup` instance.

        Parameters
        ----------
        edges : :class:`iterable` of :class:`.Edge`
            The edges in the group.

        """

        self._edge_ids = tuple(edge.get_id() for edge in edges)

    def get_edge_ids(self):
        """
        Yield the ids of edges in the group.

        Yields
        ------
        :class:`int`
            The id of an edge.

        """

        yield from self._edge_ids

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f'{self.__class__.__name__}(edge_ids={self._edge_ids})'
