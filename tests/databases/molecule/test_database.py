from tests.utilities import is_equivalent_molecule


def test_database(case_data):
    """
    Test a :class:`.MoleculeCache`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the database to test and the molecule to
        place into the database.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_database(
        database=case_data.database,
        molecule=case_data.molecule,
        key=case_data.key,
    )


def _test_database(database, molecule, key):
    """
    Test a molecule or constructed molecule database.

    Parameters
    ----------
    database : :class:`.MoleculeCache`
        The database to test.

    molecule : :class:`.Molecule`
        The molecule to put and get from the `database`.

    key : :class:`object`
        The key used to retrieve `molecule` from the database.

    Returns
    -------
    None : :class:`NoneType`

    """

    database.put(molecule)
    retrieved = database.get(key)
    is_equivalent_molecule(
        molecule1=molecule.with_canonical_atom_ordering(),
        molecule2=retrieved.with_canonical_atom_ordering(),
    )
