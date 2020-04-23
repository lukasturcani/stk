def test_get_bond(case_data):
    """
    Test :meth:`.BondInfo.get_bond`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the bond info to test and the bond it
        should be holding.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert case_data.bond_info.get_bond() is case_data.bond
