from .edge_data import EdgeData


class EdgeData_(EdgeData):
    """
    An implementation of the :class:`.EdgeData` interface.

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

    def clone(
        self,
        vertex_map=None,
        recalculate_position=False,
        add_to_vertices=True
    ):
        """
        Return a clone.

        Parameters
        ----------
        vertex_map : :class:`dict`
            If the clone should hold different :class:`.VertexData`
            instances, then a :class:`dict` should be provided, which
            maps vertex data in the current :class:`.EdgeData` to the
            vertex data instances which should be used in the clone.
            Only vertex data instances which need to be changed need
            to be present in the `vertex_map`.

        recalculate_position : :class:`bool`, optional
            Toggle if the position of the clone should be recalculated
            from the vertices it connects or if it should inherit
            the position of the original edge.

        add_to_vertices : :class:`bool`, optional
            Toggles if the clone should be added to
            :attr:`.VertexData.edges`.

        Returns
        -------
        :class:`EdgeData`
            The clone.

        """

        clone = self.__class__.__new__(self.__class__)
        clone.id = self.id
        clone.vertices = tuple(
            vertex_map.get(v, v) for v in self.vertices
        )
        clone.periodicity = np.array(self.periodicity)
        clone.custom_position = self.custom_position
        clone.lattice_constants = tuple(
            np.array(constant) for constant in self.lattice_constants
        )
        if recalculate_position:
            vertex_positions = (
                vertex.position for vertex in clone.vertices
            )
            clone.position = np.divide(
                sum(vertex_positions),
                len(clone.vertices)
            )
            self.custom_position = False
        else:
            clone.position = np.array(self.position)

        if add_to_vertices:
            for vertex in clone.vertices:
                vertex.edges.append(clone)

        return clone

    def get_edge(self):
        """
        Get an :class:`.Edge` from the data.

        Returns
        -------
        :class:`Edge`
            The edge.

        """

        return Edge(self)

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
