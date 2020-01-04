import numpy as np


class VertexData:
    """
    Holds the data used to initialize a :class:`.Vertex`.

    Attributes
    ----------
    id : :class:`int`
        The id of the vertex. Must match the index in
        :attr:`TopologyGraph.vertices`.

    position : :class:`numpy.ndarray`
        The position of the vertex.

    edges : :class:`list` of :class:`.EdgeData`
        The edges connected to the vertex.

    cell : :class:`numpy.ndarray`
        The unit cell in which the vertex is found.

    """

    def __init__(self, x, y, z):
        """
        Initialize a :class:`.VertexData` instance.

        Parameters
        ----------
        x : :class:`float`
            The x coordinate.

        y : :class:`float`
            The y coordinate.

        z : :class:`float`
            The z coordinate.

        """

        # Set by TopologyGraph.__init__.
        self.id = None
        self.position = np.array([x, y, z], dtype=np.dtype('float64'))
        self.edges = []
        self.cell = np.array([0, 0, 0])

    @classmethod
    def init_at_center(cls, *vertex_data):
        """
        Initialize at the center of other vertices.

        Parameters
        ----------
        *vertex_data : :class:`.VertexData`
            Vertices at whose center this vertex should be initialized.

        Returns
        -------
        :class:`.VertexData`
            The vertex.

        """

        center = sum(vertex.position for vertex in vertex_data)
        center /= len(vertex_data)
        return cls(*center)

    def get_vertex(self):
        """
        Get a vertex from the data.

        Returns
        -------
        :class:`.Vertex`
            A vertex initialized from the data.

        Raises
        ------
        :class:`NotImplementedError`
            This is a virtual method which needs to be implemented in
            a subclass. The implementation should return an instance of
            the matching :class:`.Vertex` subclass.

        """

        raise NotImplementedError()

    def clone(self, clear_edges=False):
        """
        Return a clone.

        Parameters
        ----------
        clear_edges : :class:`bool`, optional
            ``True`` if the clone should not be connected to any edges.

        Returns
        -------
        :class:`VertexData`
            The clone.

        """

        clone = self.__class__.__new__(self.__class__)
        for name, value in self.__dict__.items():
            if not name.startswith('_'):
                setattr(clone, name, value)

        clone.position = list(self.position)
        clone.edges = [] if clear_edges else list(self.edges)
        clone.cell = np.array(self.cell)
        return clone

    def __str__(self):
        position = self.position.tolist()
        cell_id = self.cell.tolist()
        cell = '' if cell_id == [0, 0, 0] else f', cell={cell_id}'
        return f'VertexData(id={self.id}, position={position}{cell})'

    def __repr__(self):
        return str(self)
