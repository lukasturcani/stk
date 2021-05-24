class CaseData:
    """
    A test case.

    Attributes
    ----------
    database : :class:`.MoleculeDatabase
        The database to test.

    molecule : :class:`.Molecule`
        The molecule to put and get from the :attr:`.database`.

    key : :class:`object`
        The key used to retrieve the :attr:`.molecule` from the
        database.

    """

    def __init__(self, database, molecule, key):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        database : :class:`.MoleculeDatabase`
            The database to test.

        molecule : :class:`.Molecule`
            The molecule to put and get from the `database`.

        key : :class:`object`
            The key used to retrieve the `molecule` from the database.

        """

        self.database = database
        self.molecule = molecule
        self.key = key
