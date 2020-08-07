class CaseData:
    """
    A :class:`.PeriodicInfo` test case.

    Attributes
    ----------
    periodic_info : :class:`.PeriodicInfo`
        The information to test.

    cell : :class:`tuple` of :class:`np.array` or :class:`NoneType`
        Tuple of cell lattice vectors (shape: (3,)) in Angstrom.
        `None` if topology graph is not periodic.

    lengths_and_angles : :class:`tuple` of :class:`float`
        Lengths and angles that define cell matrix,

    """

    def __init__(self, periodic_info, cell, lengths_and_angles):
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
        self.cell = cell
        self.lengths_and_angles = lengths_and_angles
