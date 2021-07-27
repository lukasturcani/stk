class CaseData:
    """
    A test case.

    Attributes
    ----------
    inchi_database : :class:`.MoleculeDatabase`
        A database to test.

    smiles_database : :class:`.MoleculeDatabase`
        A database to test.

    inchi_key_database : :class:`.MoleculeDatabase`
        A database to test.

    inchi_and_smiles_database : :class:`.MoleculeDatabase`
        The database to test collection from.

    molecules : :class:`iterable` of :class:`.Molecule`
        The molecules to put and get from the databases.

    """

    def __init__(
        self,
        inchi_database,
        smiles_database,
        inchi_key_database,
        inchi_and_smiles_database,
        molecules,
    ):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        inchi_database : :class:`.MoleculeDatabase`
            A database to test.

        smiles_database : :class:`.MoleculeDatabase`
            A database to test.

        inchi_key_database : :class:`.MoleculeDatabase`
            A database to test.

        inchi_and_smiles_database : :class:`.MoleculeDatabase`
            The database to test collection from.

        molecules : :class:`iterable` of :class:`.Molecule`
            The molecules to put and get from the databases.

        """

        self.inchi_database = inchi_database
        self.smiles_database = smiles_database
        self.inchi_key_database = inchi_key_database
        self.inchi_and_smiles_database = inchi_and_smiles_database
        self.molecules = molecules
