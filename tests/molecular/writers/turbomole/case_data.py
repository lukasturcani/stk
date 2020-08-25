class CaseData:
    """
    A :class:`.TurbomoleWriter` test case.

    Attributes
    ----------
    molecule : :class:`.Molecule`
        Molecule to test.

    writer : :class:`.TurbomoleWriter`
        The writer to test.

    string : :class:`str`
        The expected output string.

    periodic_info : :class:`.PeriodicInfo`
        Information about periodic cell. Required for testing writing
        of constructed molecules with periodic unit cells.

    """

    def __init__(self, molecule, writer, string, periodic_info):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        molecule : :class:`.Molecule`
            Molecule to test.

        writer : :class:`.TurbomoleWriter`
            The writer to test.

        string : :class:`str`
            The expected output string.

        periodic_info : :class:`.PeriodicInfo`
            Information about periodic cell.

        """

        self.molecule = molecule
        self.writer = writer
        self.string = string
        self.periodic_info = periodic_info
