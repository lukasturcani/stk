def test_get_atom1(case_data):
    """
    Test :meth:`.Bond.get_atom1`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the bond to test and the correct *atom1*.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert case_data.bond.get_atom1() is case_data.atom1
