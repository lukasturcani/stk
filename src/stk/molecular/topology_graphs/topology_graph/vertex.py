import numpy as np


class Vertex:
    """
    Represents a vertex in a :class:`.TopologyGraph`.

    """

    def __init__(self, id, position):
        """
        Initialize a :class:`.Vertex`.

        Parameters
        ----------
        id : :class:`int`
            The id of the vertex.

        position : :class:`numpy.ndarray`
            The position of the vertex.

            """

        self._id = id
        self._position = np.array(position, dtype=np.float64)

    def get_id(self):
        """
        Get the id.

        Returns
        -------
        :class:`int`
            The id.

        """

        return self._id

    def _with_scale(self, scale):
        self._position *= scale
        return self

    def with_scale(self, scale):
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

        return self.clone()._with_scale(scale)

    def clone(self):
        """
        Return a clone.

        Returns
        -------
        :class:`Vertex`
            The clone.

        """

        clone = self.__class__.__new__(self.__class__)
        clone._id = self._id
        clone._position = np.array(self._position)
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

    def _with_position(self, position):
        self._position = np.array(position, dtype=np.float64)
        return self

    def with_position(self, position):
        """
        Return a clone at a certain position.

        Parameters
        ----------
        position : :class:`numpy.ndarray`
            The desired position of the clone.

        Returns
        -------
        :class:`.Vertex`
            The clone.

        """

        return self.clone()._with_position(position)

    def get_cell(self):
        """
        Get the cell of the lattice in which the vertex is found.

        Returns
        -------
        :class:`numpy.ndarray`
            The cell of the lattice in which the vertex is found.

        """

        return np.array([0, 0, 0])

    def place_building_block(self, building_block, edges):
        """
        Place `building_block` on the :class:`.Vertex`.

        Parameters
        ----------
        building_block : :class:`.Molecule`
            The building block molecule which is to be placed on the
            vertex.

        edges : :class:`tuple` of :class:`.Edge`
            The edges to which the vertex is attached.

        Returns
        -------
        :class:`numpy.nadarray`
            The position matrix of `building_block` after being
            placed.

        """

        raise NotImplementedError()

    def map_functional_groups_to_edges(self, building_block, edges):
        """
        Map functional groups to edges.

        Each functional group in `building_block` needs to be assigned
        to an edge in `edges`.

        Parameters
        ----------
        building_block : :class:`.Molecule`
            The building block molecule which is needs to have
            functional groups assigned to edges.

        edges : :class:`tuple` of :class:`.Edge`
            The edges to which the vertex is attached.

        Returns
        -------
        :class:`dict`
            A mapping from the id of a functional group in
            `building_block` to the id of the edge in `edges` it
            is assigned to.

        """

        raise NotImplementedError()

    def __str__(self):
        position = self._position.tolist()
        return f'Vertex(id={self._id}, position={position})'

    def __repr__(self):
        return str(self)
