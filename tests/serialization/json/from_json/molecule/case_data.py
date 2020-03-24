class CaseData:
    """
    A test case.

    Attributes
    ----------
    dejsonizer : :class:`.MoleculeDejsonizer`
        The dejsonizer to test.

    json : :class:`dict`
        The JSON to test.

    position_matrix : :class:`list`
        The position matrix of the dejsonized molecule.

    molecule : :class:`.Molecule`
        The correct dejsonized molecule.

    """

    def __init__(self, dejsonizer, json, position_matrix, molecule):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        dejsonizer : :class:`.MoleculeDejsonizer`
            The dejsonizer to test.

        json : :class:`dict`
            The JSON to test.

        position_matrix : :class:`list`
            The position matrix of the dejsonized molecule.

        molecule : :class:`.Molecule`
            The correct dejsonized molecule.

        """

        self.dejsonizer = dejsonizer
        self.json = json
        self.position_matrix = position_matrix
        self.molecule = molecule
