class CaseData:
    """
    A :class:`.PeriodicInfo` test case.

    Attributes
    ----------
    periodic_info : :class:`.PeriodicInfo`
        The information to test.

    x_vector : :class:`numpy.ndarray`
        Cell lattice vector of shape (3, ) in x direction in
        Angstrom.

    y_vector : :class:`numpy.ndarray`
        Cell lattice vector of shape (3, ) in y direction in
        Angstrom.

    z_vector : :class:`numpy.ndarray`
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
        x_vector,
        y_vector,
        z_vector,
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

        cell : :class:`tuple` of :class:`np.array` or :class:`NoneType`
            Tuple of cell lattice vectors (shape: (3,)) in Angstrom.
            `None` if topology graph is not periodic.

        lengths_and_angles : :class:`tuple` of :class:`float`
            Lengths and angles that define cell matrix,

        """

        self.periodic_info = periodic_info
        self.x_vector = x_vector
        self.y_vector = y_vector
        self.z_vector = z_vector
        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
