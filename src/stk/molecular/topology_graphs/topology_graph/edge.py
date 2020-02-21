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
        Initialize a :class:`.Edge` instance.

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
            `vertex2`.

        position : :class:`numpy.ndarray`, optional
            The position of the edge, if ``None``, the midpoint of the
            the vertices is used.

        """

        if position is None:
            position = (
                vertex1.get_position() + vertex2.get_position()
            ) / 2

        self._id = id
        self._vertex1_id = vertex1.get_id()
        self._vertex2_id = vertex2.get_id()
        self._position = position
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
            The periodicity of the edge. If ``(0, 0, 0)`` the edge is
            not periodic, if ``(1, 0, -1)`` the edge is periodic going
            in the positive direction along the x axis, is not periodic
            across the y axis and is periodic in the negative direction
            along the z axis.

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

    def with_scale(self, scale):
        """
        Return a clone with a scaled position.

        Parameters
        ----------
        scale : :class:`float` or :class:`list`of :class:`float`
            The value by which the position of
            the :class:`Edge` is scaled. Can be a single number if all
            axes are scaled by the same amount or a :class:`list` of
            three numbers if each axis is scaled by a different value.

        Returns
        -------
        :class:`Edge`
            The clone.


        """

        raise NotImplementedError()

    def clone(self):
        """
        Return a clone.

        Returns
        -------
        :class:`Edge`
            The clone.

        """

        clone = self.__class__.__new__(self.__class__)
        clone._id = self._id
        clone._vertex1_id = self._vertex1_id
        clone._vertex2_id = self._vertex2_id
        clone._position = np.array(self._position)
        clone._periodicity = self._periodicity
        return clone

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
        Return the position.

        Returns
        -------
        :class:`numpy.ndarray`
            The position of the :class:`Edge`.

        """

        return np.array(self._position)
