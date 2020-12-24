class CaseData:
    """
    A :class:`.PeriodicInfo` test case.

    Attributes
    ----------
    periodic_info : :class:`.PeriodicInfo`
        The information to test.

    vector_1 : :class:`numpy.ndarray`
        Cell lattice vector of shape (3, ) in x direction in
        Angstrom.

    vector_2 : :class:`numpy.ndarray`
        Cell lattice vector of shape (3, ) in y direction in
        Angstrom.

    vector_3 : :class:`numpy.ndarray`
        Cell lattice vector of shape (3, ) in z direction in
        Angstrom.

    a : :class:`float`
        Length in a direction of cell matrix.

    b : :class:`float`
        Length in b direction of cell matrix.

    c : :class:`float`
        Length in c direction of cell matrix.

    alpha : :class:`float`
        Alpha angle of cell matrix.

    beta : :class:`float`
        Beta angle of cell matrix.

    gamma : :class:`float`
        Gamma angle of cell matrix.

    """

    def __init__(
        self,
        periodic_info,
        vector_1,
        vector_2,
        vector_3,
        a,
        b,
        c,
        alpha,
        beta,
        gamma,
    ):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        periodic_info : :class:`.PeriodicInfo`
            The information to test.

        vector_1 : :class:`numpy.ndarray`
            Cell lattice vector of shape (3, ) in x direction in
            Angstrom.

        vector_2 : :class:`numpy.ndarray`
            Cell lattice vector of shape (3, ) in y direction in
            Angstrom.

        vector_3 : :class:`numpy.ndarray`
            Cell lattice vector of shape (3, ) in z direction in
            Angstrom.

        a : :class:`float`
            Length in a direction of cell matrix.

        b : :class:`float`
            Length in b direction of cell matrix.

        c : :class:`float`
            Length in c direction of cell matrix.

        alpha : :class:`float`
            Alpha angle of cell matrix.

        beta : :class:`float`
            Beta angle of cell matrix.

        gamma : :class:`float`
            Gamma angle of cell matrix.

        """

        self.periodic_info = periodic_info
        self.vector_1 = vector_1
        self.vector_2 = vector_2
        self.vector_3 = vector_3
        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
