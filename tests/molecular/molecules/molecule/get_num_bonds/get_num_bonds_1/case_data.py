class CaseData:
    """
    A test case.

    Attributes
    ----------
    molecule : :class:`.Molecule
        The molecule to test.

    num_bonds : :class:`int`
        The correct number of bonds of :attr:`.molecule`.

    """

    def __init__(self, molecule, num_bonds):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        molecule : :class:`.Molecule`
            The molecule to test.

        num_bonds : :class:`int`
            The correct number of bonds of `molecule`.

        """

        self.molecule = molecule
        self.num_bonds = num_bonds
