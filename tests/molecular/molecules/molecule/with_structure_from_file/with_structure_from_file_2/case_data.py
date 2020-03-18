class CaseData:
    """
    A test case.

    Attributes
    ----------
    molecule : :class:`.Molecule`
        The molecule to test.

    path : :class:`str`
        The name of the structure file to use for the test.

    """

    def __init__(self, molecule, path):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        molecule : :class:`.Molecule`
            The molecule to test.

        path : :class:`str`
            The name of the structure file to use for the test.

        """

        self.molecule = molecule
        self.path = path
