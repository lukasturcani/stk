class EdgeData:
    """
    Holds data used to initialize a :class:`.Edge`.

    """

    def clone(self):
        """
        Return a clone.

        Returns
        -------
        :class:`EdgeData`
            The clone.

        """

        raise NotImplementedError()

    def with_vertices(self, vertex_map):
        """
        Return a clone holding specific vertices.

        Parameters
        ----------
        vertex_map : :class:`dict`
            Maps the id of a vertex held by the current edge to the
            :class:`.VertexData` the clone should hold instead. If
            the id of a vertex is missing, that vertex is not replaced
            in the clone.

        Returns
        -------
        :class:`.EdgeData`
            The clone.

        """

        raise NotImplementedError()

    def get_vertices(self):
        """
        Yield the vertices connected by the edge.

        Yields
        ------
        :class:`.VertexData`
            A vertex connected by the edge.

        """

        raise NotImplementedError()

    def get_periodicity(self):
        """
        Get the periodicity of the edge.

        Returns
        -------
        :class:`tuple`
            The periodicity of the edge along each dimesion.
            For example, if ``(0, 0, 0)``
            then the edge is not periodic. If, ``(1, 0, -1)`` then the
            edge is periodic across the x axis in the positive
            direction, is not periodic across the y axis and is
            periodic across the z axis in the negative direction.

        """

        raise NotImplementedError()

    def get_edge(self):
        """
        Get an :class:`.Edge` from the data.

        Returns
        -------
        :class:`Edge`
            The edge.

        """

        raise NotImplementedError()

    def get_position(self):
        """"
        Return the position of the edge.

        Returns
        -------
        :class:`numpy.ndarray`
            The position.

        """

        raise NotImplementedError()
