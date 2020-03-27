class CaseData:
    """
    A test case.

    Attributes
    ----------
    database : class:`.MoleculeValueDatabase` or \
            :class:`.ConstructedMoleculeValueDatabase`
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
        database : class:`.MoleculeValueDatabase` or \
                :class:`.ConstructedMoleculeValueDatabase`
            The database to test.

        molecule : :class:`.Molecule`
            The molecule to test.

        value : :class:`object`
            The value to put into the database.

        """

        self.database = database
        self.molecule = molecule
        self.value = value
