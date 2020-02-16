
class Vertex:
    """
    Represents a vertex in a :class:`.TopologyGraph`.

    Attributes
    ----------
    id : :class:`int`
        The id of the vertex. This should be its index in
        :attr:`TopologyGraph.vertices`.

    """

    def __init__(self, data):
        """
        Initialize a :class:`.Vertex`.

        Parameters
        ----------
        data : :class:`.VertexData`
            The vertex data.

        """

        self.id = data.id
        self._position = np.array(data.position)
        # Holds the ids of edges the Vertex is connected to.
        self._edge_ids = [edge_data.id for edge_data in data.edges]
        self._cell = np.array(data.cell)
        # This holds the ConstructedMolecule that the vertex is used
        # to construct.
        self._mol = None

    def apply_scale(self, scale):
        """
        Scale the position by `scale`.

        Parameters
        ----------
        scale : :class:`float` or :class:`list`of :class:`float`
            The value by which the position of the :class:`Vertex` is
            scaled. Can be a single number if all axes are scaled by
            the same amount or a :class:`list` of three numbers if
            each axis is scaled by a different value.

        Returns
        -------
        :class:`Vertex`
            The vertex is returned.

        """

        self._position *= scale
        return self

    def clone(self, clear_edges=False):
        """
        Return a clone.

        Parameters
        ----------
        clear_edges : :class:`bool`, optional
            ``True`` if the clone should not be connected to any edges.

        Returns
        -------
        :class:`Vertex`
            The clone.

        """

        clone = self.__class__.__new__(self.__class__)
        for name, value in self.__dict__.items():
            if not name.startswith('_'):
                setattr(clone, name, value)
        clone._position = np.array(self._position)
        clone._cell = np.array(self._cell)
        clone._edge_ids = [] if clear_edges else list(self._edge_ids)
        clone._mol = self._mol
        return clone

    def get_position(self):
        """
        Return the position.

        Returns
        -------
        :class:`numpy.ndarray`
            The position of the :class:`Vertex`.

        """

        return np.array(self._position)

    def get_num_edges(self):
        """
        Return the number of connceted edge.

        Returns
        -------
        :class:`int`
            The number of connected edges.

        """

        return len(self._edge_ids)

    def get_edge_ids(self):
        """
        Yield the ids of connected edges.

        Yields
        ------
        :class:`int`
            The :class:`~.Edge.id` of a connected edge.

        """

        yield from self._edge_ids

    def get_cell(self):
        """
        Get the cell of the lattice in which the vertex is found.

        Returns
        -------
        :class:`numpy.ndarray`
            The cell of the lattice in which the vertex is found.

        """

        return np.array(self._cell)

    def set_constructed_molecule(self, mol):
        """
        Set the :class:`.ConstructedMolecule` being constructed.

        Parameters
        ----------
        mol : :class:`.ConstructedMolecule`
            The molecule being constructed.

        Returns
        -------
        :class:`.Vertex`
            The vertex.

        """

        self._mol = mol
        return self

    def place_building_block(self, building_block, vertices, edges):
        """
        Place `building_block` on the :class:`.Vertex`.

        Parameters
        ----------
        building_block : :class:`.Molecule`
            The building block molecule which is to be placed on the
            vertex.

        vertices : :class:`tuple` of :class:`.Vertex`
            All vertices in the topology graph. The index of each
            vertex must match its :class:`~.Vertex.id`.

        edges : :class:`tuple` of :class:`.Edge`
            All edges in the topology graph. The index of each
            edge must match its :class:`~.Edge.id`.

        Returns
        -------
        :class:`numpy.nadarray`
            The position matrix of `building_block` after being
            placed.

        Raises
        ------
        :class:`NotImplementedError`
            This is a virtual method, it needs to be implemented in a
            subclass.

        """

        raise NotImplementedError()

    def assign_func_groups_to_edges(
        self,
        building_block,
        vertices,
        edges
    ):
        """
        Assign functional groups to edges.

        Each :class:`.FunctionalGroup` of the `building_block` needs
        to be associated with one of the :class:`.Edge` instances in
        :attr:`edges`.

        Parameters
        ----------
        building_block : :class:`.Molecule`
            The building block molecule which is needs to have
            functional groups assigned to edges.

        vertices : :class:`tuple` of :class:`.Vertex`
            All vertices in the topology graph. The index of each
            vertex must match its :class:`~.Vertex.id`.

        edges : :class:`tuple` of :class:`.Edge`
            All edges in the topology graph. The index of each
            edge must match its :class:`~.Edge.id`.

        Returns
        -------
        :class:`dict`
            A mapping from the id of a functional group in
            `building_block` to the id of the edge in :attr:`edges` it
            is assigned to.

        Raises
        ------
        :class:`NotImplementedError`
            This is a virtual method, it needs to be implemented in a
            subclass.

        """

        raise NotImplementedError()

    def after_assign_func_groups_to_edges(
        self,
        building_block,
        func_groups,
        vertices,
        edges
    ):
        """
        Perform operations after functional groups have been assigned.

        This method is always executed serially. It is often useful
        when data needs to be transferred between vertices, which
        have been processed independently, in parallel.

        It does nothing by default, but should be overridden when
        necessary.

        Parameters
        ----------
        building_block : :class:`.Molecule`
            The building block molecule which is needs to have
            functional groups assigned to edges.

        func_groups : :class:`tuple` of :class:`.FunctionalGroup`
            The functional group clones added to the constructed
            molecule.

        vertices : :class:`tuple` of :class:`.Vertex`
            All vertices in the topology graph. The index of each
            vertex must match its :class:`~.Vertex.id`.

        edges : :class:`tuple` of :class:`.Edge`
            All edges in the topology graph. The index of each
            edge must match its :class:`~.Edge.id`.

        Returns
        -------
        None : :class:`NoneType`

        """

        return

    def _get_edge_centroid(self, centroid_edges, vertices):
        """
        Return the centroid of `centroid_edges`.

        Parameters
        ----------
        centroid_edges : :class:`iterable` of :class:`.Edge`
            The edges which are used to calculate the centroid.

        vertices : :class:`tuple` of :class:`.Vertex`
            All the vertices in the topology graph. Index of each
            vertex must be equal to its :class:`~.Vertex.id`.

        Returns
        -------
        :class:`numpy.ndarray`
            The centroid of the edges.

        """

        edge_positions = []
        for i, edge in enumerate(centroid_edges, 1):
            edge_positions.append(edge.get_position(self, vertices))
        return np.sum(edge_positions, axis=0) / i

    def _get_edge_plane_normal(self, reference, plane_edges, vertices):
        """
        Get the normal to the plane on which `edges` lie.

        Parameters
        ----------
        reference : :class:`numpy.ndarray`
            A reference direction vector. The direction of the returned
            normal is set such that its angle with with `reference`
            is always acute.

        plane_edges : :class:`iterable` of :class:`.Edge`
            The edges which are used to calculate the plane.
            If there are more than three, a plane of best fit across
            `edges` is returned.

        vertices : :class:`tuple` of :class:`.Vertex`
            All the vertices in the topology graph. Index of each
            vertex must be equal to its :class:`~.Vertex.id`.

        Returns
        -------
        :class:`numpy.ndarray`
            A unit vector which describes the normal to the plane of
            the edges.

        """

        edge_positions = []
        for i, edge in enumerate(plane_edges, 1):
            edge_positions.append(edge.get_position(self, vertices))
        edge_positions = np.array(edge_positions)

        centroid = np.sum(edge_positions, axis=0) / i
        normal = np.linalg.svd(edge_positions - centroid)[-1][2, :]

        if vector_angle(normal, reference) > np.pi/2:
            normal *= -1
        return normal

    def _get_molecule_centroid(self, atom_ids=None):
        """
        Get the centroid of the molecule being constructed.

        During construction :meth:`.Molecule.get_centroid` cannot be
        used, because the molecule is not fully constructed yet. This
        method acts as its replacement during construction.

        Parameters
        ----------
        atom_ids : :class:`iterable` of :class:`int`, optional
            The ids of atoms which are used to calculate the
            centroid. If ``None``, then all atoms will be used.

        Returns
        -------
        :class:`numpy.ndarray`
            The centroid of atoms specified by `atom_ids`.

        """

        if atom_ids is None:
            atom_ids = range(len(self._mol.atoms))
        elif not isinstance(atom_ids, (list, tuple)):
            atom_ids = list(atom_ids)

        return np.divide(
            np.sum(
                np.array(self._mol._position_matrix)[atom_ids, :],
                axis=0
            ),
            len(atom_ids)
        )

    def __str__(self):
        position = self._position.tolist()
        cell_id = self._cell.tolist()
        cell = '' if cell_id == [0, 0, 0] else f', cell={cell_id}'
        return f'Vertex(id={self.id}, position={position}{cell})'

    def __repr__(self):
        return str(self)
