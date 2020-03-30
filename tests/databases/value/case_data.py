class CaseData:
    """
    A test case.

    Attributes
    ----------
    database : class:`.ValueDatabase`
        The database to test.

    molecule : :class:`.Molecule`
        The molecule to test.

    value : :class:`object`
        The value to put into the database.

    """

    def __init__(self, database, molecule, value):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        database : class:`.ValueDatabase`
            The database to test.

        molecule : :class:`.Molecule`
            The molecule to test.

        value : :class:`object`
            The value to put into the database.

        """

        self.database = database
        self.molecule = molecule
        self.value = value
