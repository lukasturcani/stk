class CaseData:
    """
    A test case.

    Attributes
    ----------
    start : :class:`numpy.ndarray`
        The starting position of the tested edge.

    scale : :class:`float` or :class:`numpy.ndarray`
        The scale applied to the tested edge.

    target : :class:`numpy.ndarray`
        The expected position of the new edge.

    """

    def __init__(self, start, scale, target):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        start : :class:`numpy.ndarray`
            The starting position of the tested edge.

        scale : :class:`float` or :class:`numpy.ndarray`
            The scale applied to the tested edge.

        target : :class:`numpy.ndarray`
            The expected position of the new edge.

        """

        self.start = start
        self.scale = scale
        self.target = target
