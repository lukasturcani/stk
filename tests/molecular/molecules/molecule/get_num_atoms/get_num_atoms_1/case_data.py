class CaseData:
    """
    A test case.

    Attributes
    ----------
    molecule : :class:`.Molecule`
        The molecule to test.

    num_atoms : :class:`int`
        The correct number of atoms :attr:`.molecule` should have.

    """

    def __init__(self, molecule, num_atoms):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        molecule : :class:`.Molecule`
            The molecule to test.

        num_atoms : :class:`int`
            The correct number of atoms `molecule` should have.

        """

        self.molecule = molecule
        self.num_atoms = num_atoms
