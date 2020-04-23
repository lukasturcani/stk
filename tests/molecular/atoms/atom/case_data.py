class CaseData:
    """
    An :class:`.Atom` test case.

    Attributes
    ----------
    atom : :class:`.Atom`
        The atom being tested.

    id : :class:`int`
        The correct id.

    charge : :class:`int`
        The correct charge.

    atomic_number : :class:`int`
        The correct atomic number.

    """

    def __init__(self, atom, id, charge, atomic_number):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        atom : :class:`.Atom`
            The atom being tested.

        id : :class:`int`
            The correct id.

        charge : :class:`int`
            The correct charge.

        atomic_number : :class:`int`
            The correct atomic number.

        """

        self.atom = atom
        self.id = id
        self.charge = charge
        self.atomic_number = atomic_number
