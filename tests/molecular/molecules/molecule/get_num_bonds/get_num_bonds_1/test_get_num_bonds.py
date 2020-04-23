def test_get_num_bonds(case_data):
    """
    Test :meth:`.Molecule.get_num_bonds`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the molecule to test and the correct number
        of bonds.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert case_data.molecule.get_num_bonds() == case_data.num_bonds
