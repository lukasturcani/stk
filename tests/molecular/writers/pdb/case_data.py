import stk


class CaseData:
    """
    A :class:`.PdbWriter` test case.

    Attributes:
        molecule:
            Molecule to test.

        writer:
            The writer to test.

        string:
            The expected output string.

        periodic_info:
            Information about periodic cell. Required for testing writing
            of constructed molecules with periodic unit cells.

    """

    def __init__(
        self,
        molecule: stk.Molecule,
        writer: stk.PdbWriter,
        string: str,
        periodic_info: stk.PeriodicInfo | None,
    ) -> None:
        """
        Initialize a :class:`.CaseData` instance.

        Parameters:
            molecule:
                Molecule to test.

            writer:
                The writer to test.

            string:
                The expected output string.

            periodic_info:
                Information about periodic cell.

        """

        self.molecule = molecule
        self.writer = writer
        self.string = string
        self.periodic_info = periodic_info
