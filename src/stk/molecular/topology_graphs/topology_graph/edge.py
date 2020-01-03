
class Edge:
    """
    Represents an edge in a topology graph.

    Note that some methods of this class will behave differently
    before and after :meth:`finalize` is called. Methods will switch
    from returning :class:`.Vertex` objects to returning :class:`int`
    objects.

    Attributes
    ----------
    id : :class:`int`
        The id of the edge. Matches the index of the edge in
        :attr:`.TopologyGraph.edges`.

    """

    def __init__(self, data):
        """
        Initialize an :class:`Edge`.

        Parameters
        ----------
        data : :class:`.EdgeData`
            The edge data.

        """

        self._vertex_ids = [vertex.id for vertex in data.vertices]
        self.id = data.id
        self._periodicity = np.array(data.periodicity)
        # The FunctionalGroup instances which the edge connects.
        # These will belong to the molecules placed on the vertices
        # connected by the edge.
        self._func_groups = []

        self._custom_position = data.custom_position
        self._position = np.array(data.position)
        self._lattice_constants = tuple(
            np.array(constant) for constant in data.lattice_constants
        )

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

        return np.array(self._periodicity)

    def is_periodic(self):
        """
        Return ``True`` if periodic.

        Returns
        -------
        :class:`bool`
            ``True`` if periodic.

        """

        return any(i != 0 for i in self._periodicity)

    def apply_scale(self, scale):
        """
        Scale the position by `scale`.

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
            The edge is returned.


        """

        self._position *= scale
        self._lattice_constants = tuple(
            scale*constant for constant in self._lattice_constants
        )
        return self

    def clone(self):
        """
        Return a clone.

        Returns
        -------
        :class:`Edge`
            The clone.

        """

        clone = self.__class__.__new__(self.__class__)
        clone.id = self.id
        clone._func_groups = list(self._func_groups)
        clone._custom_position = self._custom_position
        clone._periodicity = np.array(self._periodicity)
        clone._lattice_constants = tuple(
            np.array(constant) for constant in self._lattice_constants
        )
        clone._vertex_ids = tuple(self._vertex_ids)
        clone._position = np.array(self._position)
        return clone

    def get_func_groups(self):
        """
        Get the functional groups connected by this edge.

        Returns
        -------
        :class:`tuple` of :class:`.FunctionalGroup`
            The functional groups connected by the edge.

        """

        return tuple(self._func_groups)

    def get_vertex_ids(self):
        """
        Get the connected vertices.

        Yields
        ------
        :class:`int`
            The id of a connected vertex.

        """

        yield from self._vertex_ids

    def assign_func_group(self, func_group):
        """
        Assign `func_group` to be connected by this edge.

        Parameters
        ----------
        func_group : :class:`.FunctionalGroup`
            The functional group to be assigned to the edge.

        Returns
        -------
        :class:`Edge`
            The edge is returned.

        """

        self._func_groups.append(func_group)

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

        if reference is None or not self.is_periodic():
            return np.array(self._position)

        other = vertices[
            next(v for v in self._vertex_ids if v != reference.id)
        ]
        direction = (
            1 if reference is vertices[self._vertex_ids[0]] else -1
        )
        end_cell = reference.get_cell() + direction*self._periodicity
        cell_shift = end_cell - other.get_cell()
        shift = 0
        for dim, constant in zip(cell_shift, self._lattice_constants):
            shift += dim*constant
        return (
            (other.get_position()+shift+reference.get_position()) / 2
        )

    def __str__(self):
        return repr(self)

    def __repr__(self):
        vertices = ', '.join(str(id_) for id_ in self._vertex_ids)
        if self._custom_position:
            position = f', position={self._position!r}'
        else:
            position = ''

        if any(i != 0 for i in self._periodicity):
            periodicity = f', periodicity={tuple(self._periodicity)!r}'
        else:
            periodicity = ''

        return f'Edge({vertices}{position}{periodicity})'
