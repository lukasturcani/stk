import numpy as np


class EdgeData:
    """
    Holds data used to initialize a :class:`.Edge`.

    """

    def __init__(
        self,
        *vertex_data,
        position=None,
        periodicity=None,
        lattice_constants=None
    ):
        """
        Initialize an :class:`EdgeData` instance.

        Parameters
        ----------
        *vertex_data : :class:`.VertexData`
            The vertices which the :class:`Edge` connects.

        position : :class:`numpy.ndarray`, optional
            The position of the edge. If ``None``, the centroid
            of `vertex_data` is used.

        periodicity : :class:`tuple` of :class:`int`, optional
            The periodicity of the edge. For example, if ``(0, 0, 0)``
            then the edge is not periodic. If, ``(1, 0, -1)`` then the
            edge is periodic across the x axis in the positive
            direction, is not periodic across the y axis and is
            periodic across the z axis in the negative direction. If
            ``None`` then the edge is not periodic.

        lattice_constants : :class:`iterable`, optional
            If the edge is periodic, the a, b and c lattice
            constants should be provided as vectors in Cartesian
            coordinates.

        """

        if periodicity is None:
            periodicity = [0, 0, 0]
        if lattice_constants is None:
            lattice_constants = ([0, 0, 0] for i in range(3))

        # This is set by TopologyGraph.__init__.
        self.id = None
        self.vertices = vertex_data
        self.periodicity = np.array(periodicity)
        self.custom_position = position is not None
        self.position = position
        self.lattice_constants = tuple(
            np.array(constant) for constant in lattice_constants
        )

        _position = 0
        for i, vertex in enumerate(vertex_data, 1):
            vertex.edges.append(self)

            if not self.custom_position:
                _position += vertex.position

        if not self.custom_position:
            self.position = _position / i

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

    def __str__(self):
        return repr(self)

    def __repr__(self):
        vertices = ', '.join(
            str(vertex.id) for vertex in self.vertices
        )
        id_ = '' if self.id is None else f', id={self.id}'
        if self.custom_position:
            position = f', position={self.position!r}'
        else:
            position = ''

        if any(i != 0 for i in self.periodicity):
            periodicity = f', periodicity={tuple(self.periodicity)!r}'
        else:
            periodicity = ''

        return f'EdgeData({vertices}{id_}{position}{periodicity})'
