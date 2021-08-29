from __future__ import annotations

import stk


class CaseData:
    """
    A test case.

    Attributes:

        database: The database to test.

        expected_molecules: The expected molecules to get from the
            databases using their smiles as the key.

    """

    def __init__(
        self,
        database: stk.MoleculeDatabase,
        expected_molecules: dict[str, stk.BuildingBlock],
    ):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters:

            database: The database to test.

            expected_molecules: The expected molecules to get from the
                databases using their smiles as the key.

        """

        self.database = database
        self.expected_molecules = dict(expected_molecules)
