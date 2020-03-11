"""
Edge
====

"""

import numpy as np


class Edge:
    """
    Represents an edge in a topology graph.

    """

    def __init__(
        self,
        id,
        vertex1,
        vertex2,
        periodicity=(0, 0, 0),
        position=None,
    ):
        """
        Initialize an :class:`.Edge` instance.

        Parameters
        ----------
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

        if position is None:
            self._position = (
                vertex1.get_position() + vertex2.get_position()
            ) / 2
        else:
            self._position = np.array(position, dtype=np.float64)

        self._id = id
        self._vertex1_id = vertex1.get_id()
        self._vertex2_id = vertex2.get_id()
        self._periodicity = periodicity

    def get_id(self):
        """
        Get the id of the edge.

        Returns
        -------
        :class:`int`
            The id.

        """

        return self._id

    def get_periodicity(self):
        """
        Get the periodicity of the edge.

        Returns
        -------
        :class:`tuple` of :class:`int`
            The periodicity of the edge. For example, if ``(0, 0, 0)``
            the edge is not periodic, if ``(1, 0, -1)`` the edge is
            periodic going in the positive direction along the x axis,
            is not periodic across the y axis and is periodic in the
            negative direction along the z axis.

        """

        return self._periodicity

    def is_periodic(self):
        """
        Return ``True`` if periodic.

        Returns
        -------
        :class:`bool`
            ``True`` if periodic.

        """

        return any(i != 0 for i in self._periodicity)

    def _with_scale(self, scale):
        """
        Modify the edge.

        """

        self._position *= scale
        return self

    def with_scale(self, scale):
        """
        Return a clone with a scaled position.

        Parameters
        ----------
        scale : :class:`float` or :class:`tuple` of :class:`float`
            The value by which the position of
            the :class:`.Edge` is scaled. Can be a single number if all
            axes are scaled by the same amount or a :class:`tuple` of
            three numbers if each axis is scaled by a different value.

        Returns
        -------
        :class:`.Edge`
            The clone. Has the same type as the original edge.


        """

        return self.clone()._with_scale(scale)

    def clone(self):
        """
        Return a clone.

        Returns
        -------
        :class:`.Edge`
            The clone. Has the same type as the original edge.

        """

        clone = self.__class__.__new__(self.__class__)
        clone._id = self._id
        clone._vertex1_id = self._vertex1_id
        clone._vertex2_id = self._vertex2_id
        clone._position = np.array(self._position)
        clone._periodicity = self._periodicity
        return clone

    def get_vertex_ids(self):
        """
        Yield the ids of the vertices connected by the edge.

        The id of the first vertex is yielded first, followed by the id
        of the second vertex.

        Yields
        ------
        :class:`int`
            The id of a :class:`.Vertex`.

        """

        yield self._vertex1_id
        yield self._vertex2_id

    def get_vertex1_id(self):
        """
        Get the id of the first vertex.

        Returns
        -------
        :class:`int`
            The id of the first vertex.

        """

        return self._vertex1_id

    def get_vertex2_id(self):
        """
        Get the id of the second vertex.

        Returns
        -------
        :class:`int`
            The id of the second vertex.

        """

        return self._vertex2_id

    def get_position(self):
        """
        Get the position.

        Returns
        -------
        :class:`numpy.ndarray`
            The position of the :class:`Edge`.

        """

        return np.array(self._position)

    def _with_position(self, position):
        """
        Modify the edge.

        """

        self._position = np.array(position, dtype=np.float64)
        return self

    def with_position(self, position):
        """
        Return a clone at `position`.

        Parameters
        ----------
        position : :class:`numpy.ndarray`
            The desired position of the clone.

        Returns
        -------
        :class:`.Edge`
            The clone. Has the same type as the original edge.

        """

        return self.clone()._with_position(position)

    def __repr__(self):
        return str(self)

    def __str__(self):
        return (
            f'Edge({self._id}, {self._vertex1_id}, {self._vertex2_id})'
        )
