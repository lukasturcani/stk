class Edge:
    """
    Represents an edge in a topology graph.

    """

    def get_periodicity(self):
        """
        Get the periodicity of the edge.

        Returns
        -------
        :class:`numpy.ndarray`
            The periodicity of the edge. If ``[0, 0, 0]`` the edge is
            not periodic, if ``[1, 0, -1]`` the edge is periodic going
            in the positive direction along the x axis, is not periodic
            across the y axis and is periodic in the negative direction
            along the z axis.

        """

        raise NotImplementedError()

    def is_periodic(self):
        """
        Return ``True`` if periodic.

        Returns
        -------
        :class:`bool`
            ``True`` if periodic.

        """

        return any(i != 0 for i in self.get_periodicity())

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

        raise NotImplementedError()

    def get_vertex_ids(self):
        """
        Get the ids of connected vertices.

        Yields
        ------
        :class:`int`
            The id of a connected vertex.

        """

        raise NotImplementedError()

    def get_position(self, reference=None, vertices=None):
        """
        Return the position.

        Parameters
        ----------
        reference : :class:`.Vertex`, optional
            If the edge is periodic, the position returned will
            depend on which vertex the edge position is calculated
            relative to.

        vertices : :class:`tuple` of :class:`.Vertex`, optional
            All the vertices in the topology graph. Index of each
            vertex must be equal to its :class:`~.Vertex.id`. Only
            needs to be supplied if `reference` is supplied.

        Returns
        -------
        :class:`numpy.ndarray`
            The position of the :class:`Edge`.

        """

        raise NotImplementedError()
