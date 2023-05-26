from __future__ import annotations

import stk

from tests.utilities import is_equivalent_constructed_molecule

from .case_data import CaseData


def test_get_all(case_data: CaseData) -> None:
    """
    Test iteration over all entries.

    Parameters:

        case_data: A test case. Holds the databases to test and the
            molecules to place into the databases.

    """

    _test_get_all(
        database=case_data.database,
        expected_molecules=case_data.expected_molecules,
    )


def _test_get_all(
    database: stk.ConstructedMoleculeDatabase,
    expected_molecules: dict[str, stk.ConstructedMolecule],
) -> None:
    """
    Test iteration over all entries.

    Parameters:

        database: A database to test.

        expected_molecules: The expected molecules to get from the
            databases using their smiles as the key.

    """

    smiles = stk.Smiles()
    for i, retrieved in enumerate(database.get_all()):
        expected = expected_molecules[smiles.get_key(retrieved)]
        is_equivalent_constructed_molecule(
            constructed_molecule1=(expected.with_canonical_atom_ordering()),
            constructed_molecule2=(retrieved.with_canonical_atom_ordering()),
        )

    assert i + 1 == len(expected_molecules)
