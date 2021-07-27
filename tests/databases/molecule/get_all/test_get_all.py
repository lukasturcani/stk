import stk

from tests.utilities import is_equivalent_molecule


def test_get_all(case_data):
    """
    Test iteration over all entries.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the databases to test and the molecules to
        place into the databases.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_get_all(
        inchi_database=case_data.inchi_database,
        smiles_database=case_data.smiles_database,
        inchi_key_database=case_data.inchi_key_database,
        inchi_and_smiles_database=case_data.inchi_and_smiles_database,
        molecules=case_data.molecules,
    )


def _test_get_all(
    inchi_database,
    smiles_database,
    inchi_key_database,
    inchi_and_smiles_database,
    molecules,
):
    """
    Test iteration over all entries.

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

    Returns
    -------
    None : :class:`NoneType`

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
