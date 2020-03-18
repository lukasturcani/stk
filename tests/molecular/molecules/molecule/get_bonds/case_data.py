class CaseData:
    """
    A test case.

    Attributes
    ----------
    molecule : :class:`.Molecule`
        The molecule to test.

    bonds : :class:`tuple` of :class:`.Bond`
        The correct bonds of :attr:`.molecule`.

    """

    def __init__(self, molecule, bonds):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        molecule : :class:`.Molecule`
            The molecule to test.

        bonds : :class:`tuple` of :class:`.Bond`
            The correct bonds of `molecule`.

        """

        self.molecule = molecule
        self.bonds = bonds
