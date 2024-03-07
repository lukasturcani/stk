import stk


class CaseData:
    """
    A :class:`.XyzWriter` test case.

    Attributes:
        molecule:
            Molecule to test.

        writer:
            The writer to test.

        string:
            The expected output string.

    """

    def __init__(
        self,
        molecule: stk.Molecule,
        writer: stk.XyzWriter,
        string: str,
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

        """

        self.molecule = molecule
        self.writer = writer
        self.string = string
