class CaseData:
    """
    A test case.

    Attributes
    ----------
    vertex : :class:`.Vertex`
        The vertex to test.

    id : :class:`int`
        The correct id of :attr:`.vertex`.

    position : :class:`numpy.ndarray`
        The correct position of :attr:`.vertex`.

    cell : :class:`numpy.ndarray`
        The correct position of :attr:`.vertex`.

    """

    def __init__(self, vertex, id, position, cell):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        vertex : :class:`.Vertex`
            The vertex to test.

        id : :class:`int`
            The correct id of :attr:`.vertex`.

        position : :class:`numpy.ndarray`
            The correct position of :attr:`.vertex`.

        cell : :class:`numpy.ndarray`
            The correct position of :attr:`.vertex`.

        """

        self.vertex = vertex
        self.id = id
        self.position = position
        self.cell = cell
