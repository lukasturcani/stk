class CaseData:
    """
    A test case.

    Attributes
    ----------
    edge : :class:`.Edge`
        The edge to test.

    id : :class:`int`
        The correct id of the edge.

    vertex1_id : :class:`int`
        The correct id of the first vertex held by the edge.

    vertex2_id : :class:`int`
        The correct id of the second vertex held by the edge.

    periodicity : :class:`tuple` of :class:`int`
        The correct periodicity of :attr:`.edge`.

    is_periodic : :class:`bool`
        The truth about :attr:`.edge` being periodic.

    """

    def __init__(
        self,
        edge,
        id,
        vertex1_id,
        vertex2_id,
        periodicity,
        is_periodic,
    ):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        edge : :class:`.Edge`
            The edge to test.

        id : :class:`int`
            The correct id of the edge.

        vertex1_id : :class:`int`
            The correct id of the first vertex held by the edge.

        vertex2_id : :class:`int`
            The correct id of the second vertex held by the edge.

        periodicity : :class:`tuple` of :class:`int`
            The correct periodicity of `edge`.

        is_periodic : :class:`bool`
            The truth about `edge` being periodic.

        """

        self.edge = edge
        self.id = id
        self.vertex1_id = vertex1_id
        self.vertex2_id = vertex2_id
        self.periodicity = periodicity
        self.is_periodic = is_periodic
