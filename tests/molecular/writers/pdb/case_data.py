class CaseData:
    """
    A :class:`.PdbWriter` test case.

    Attributes
    ----------
    molecule : :class:`.Molecule`
        Molecule to test.

    writer : :class:`.PdbWriter`
        The writer to test.

    string : :class:`str`
        The expected output string.

    periodic_cell : :class:`tuple` of :class:`np.array`
        Tuple of cell lattice vectors (shape: (3,)) in Angstrom.
        Required for testing writing of constructed molecules with
        periodic unit cells.

    """

    def __init__(self, molecule, writer, string, periodic_cell):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        molecule : :class:`.Molecule`
            Molecule to test.

        writer : :class:`.PdbWriter`
            The writer to test.

        string : :class:`str`
            The expected output string.

        periodic_cell : :class:`tuple` of :class:`np.array`
            Tuple of cell lattice vectors (shape: (3,)) in Angstrom.

        """

        self.molecule = molecule
        self.writer = writer
        self.string = string
        self.periodic_cell = periodic_cell
