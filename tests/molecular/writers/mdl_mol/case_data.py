class CaseData:
    """
    A :class:`.MolWriter` test case.

    Attributes
    ----------
    molecule : :class:`.Molecule`
        Molecule to test.

    writer : :class:`.MolWriter`
        The writer to test.

    string : :class:`str`
        The expected output string.

    """

    def __init__(self, molecule, writer, string):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        molecule : :class:`.Molecule`
            Molecule to test.

        writer : :class:`.MolWriter`
            The writer to test.

        string : :class:`str`
            The expected output string.

        """

        self.molecule = molecule
        self.writer = writer
        self.string = string
