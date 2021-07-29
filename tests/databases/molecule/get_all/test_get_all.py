from typing import Iterable
import stk

from tests.utilities import is_equivalent_molecule
from .case_data import CaseData


def test_get_all(case_data: CaseData) -> None:
    """
    Test iteration over all entries.

    Parameters:

        case_data: A test case. Holds the databases to test and the
            molecules to place into the databases.

    Returns:

        None.

    """

    _test_get_all(
        inchi_database=case_data.inchi_database,
        smiles_database=case_data.smiles_database,
        inchi_key_database=case_data.inchi_key_database,
        inchi_and_smiles_database=case_data.inchi_and_smiles_database,
        molecules=case_data.molecules,
    )


def _test_get_all(
    inchi_database: stk.MoleculeDatabase,
    smiles_database: stk.MoleculeDatabase,
    inchi_key_database: stk.MoleculeDatabase,
    inchi_and_smiles_database: stk.MoleculeDatabase,
    molecules: Iterable[stk.BuildingBlock],
) -> None:
    """
    Test iteration over all entries.

    Parameters:

        inchi_database: A database to test.

        smiles_database: A database to test.

        inchi_key_database: A database to test.

        inchi_and_smiles_database The database to test.

        molecules: The molecules to put and get from the databases.

    Returns:

        None.

    """

    smiles = stk.Smiles()

    inchi_molecules = molecules[:2]
    smiles_molecules = molecules[2:4]
    inchi_and_smiles_molecules = molecules[4:]

    for molecule in inchi_molecules:
        inchi_database.put(molecule)

    for molecule in smiles_molecules:
        smiles_database.put(molecule)

    for molecule in inchi_and_smiles_molecules:
        inchi_and_smiles_database.put(molecule)

    smiles_to_molecule = {
        smiles.get_key(molecule): molecule
        for molecule in molecules
    }

    # Use an InChIKey database for get_all because none of the
    # molecules use this key, but this should not matter.
    for i, retrieved in enumerate(inchi_key_database.get_all()):
        expected = smiles_to_molecule[smiles.get_key(retrieved)]
        is_equivalent_molecule(
            molecule1=expected.with_canonical_atom_ordering(),
            molecule2=retrieved.with_canonical_atom_ordering(),
        )

    assert i+1 == len(molecules)
