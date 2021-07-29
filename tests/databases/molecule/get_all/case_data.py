from typing import Iterable
import stk


class CaseData:
    """
    A test case.

    Attributes:

        inchi_database: A database to test.

        smiles_database: A database to test.

        inchi_key_database: A database to test.

        inchi_and_smiles_database: The database to test.

        molecules: The molecules to put and get from the databases.

    """

    def __init__(
        self,
        inchi_database: stk.MoleculeDatabase,
        smiles_database: stk.MoleculeDatabase,
        inchi_key_database: stk.MoleculeDatabase,
        inchi_and_smiles_database: stk.MoleculeDatabase,
        molecules: Iterable[stk.BuildingBlock],
    ):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters:

            inchi_database: A database to test.

            smiles_database: A database to test.

            inchi_key_database: A database to test.

            inchi_and_smiles_database The database to test.

            molecules: The molecules to put and get from the databases.

        """

        self.inchi_database = inchi_database
        self.smiles_database = smiles_database
        self.inchi_key_database = inchi_key_database
        self.inchi_and_smiles_database = inchi_and_smiles_database
        self.molecules = molecules