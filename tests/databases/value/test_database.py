def test_database(case_data):
    """
    Test a database.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the database to test and the value to put
        into it.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_database(
        database=case_data.database,
        molecule=case_data.molecule,
        value=case_data.value,
    )


def _test_database(database, molecule, value):
    """
    Test a database.

    Parameters
    ----------
    database : class:`.ValueDatabase`
        The database to test.

    molecule : :class:`.Molecule`
        The molecule to test.

    value : :class:`object`
        The value to put into the database.

    Returns
    -------
    None : :class:`NoneType`

    """

    database.put(molecule, value)
    assert database.get(molecule) == value
