"""
COF Edge
========

"""

from ..topology_graph import Edge


class CofEdge(Edge):
    """
    A :class:`.Cof` edge.

    """

    def __init__(
        self,
        parent_id,
        id,
        vertex1,
        vertex2,
        periodicity=(0, 0, 0),
        position=None,
    ):
        """
        Initialize a :class:`.CofEdge`.

        Parameters
        ----------
        parent_id : :class:`int`
            The id of the edge in the (1, 1, 1) unit cell, with which
            this one is periodic.

        id : :class:`int`
            The id of the edge.

        vertex1 : :class:`.Vertex`
            The first vertex the edge is connected to.

        vertex2 : :class:`.Vertex`
            The second vertex the edge is connected to.

        periodicity : :class:`tuple` of :class:`int`, optional
            The periodicity of the edge, when going from `vertex1` to
            `vertex2`. For example, if ``(0, 0, 0)`` the edge is not
            periodic, if ``(1, 0, -1)`` the edge is periodic going in
            the positive direction along the x axis, is not periodic
            across the y axis and is periodic in the negative direction
            along the z axis.

        position : :class:`numpy.ndarray`, optional
            The position of the edge, if ``None``, the midpoint of the
            vertices is used.

        """

        super().__init__(id, vertex1, vertex2, periodicity, position)
        self._parent_id = parent_id

    def get_parent_id(self):
        """
        Get the id of the *parent* edge.

        The *parent* edge is the edge in the (1, 1, 1) unit cell, with
        which this one is periodic.

        Returns
        -------
        :class:`int`
            The id of the parent edge.

        """

        return self._parent_id

    def clone(self):
        obj = super().clone()
        obj._parent_id = self._parent_id
        return obj
