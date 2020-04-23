class CaseData:
    """
    A test case.

    Attributes
    ----------
    dejsonizer : :class:`.MoleculeDejsonizer`
        The dejsonizer to test.

    json : :class:`dict`
        The JSON to test.

    molecule : :class:`.Molecule`
        The correct dejsonized molecule.

    """

    def __init__(self, dejsonizer, json, molecule):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        dejsonizer : :class:`.MoleculeDejsonizer`
            The dejsonizer to test.

        json : :class:`dict`
            The JSON to test.

        molecule : :class:`.Molecule`
            The correct dejsonized molecule.

        """

        self.dejsonizer = dejsonizer
        self.json = json
        self.molecule = molecule
