def test_is_dative(case_data):
    """
    Test :meth:`.Bond.is_dative`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the bond to test and the correct
        dativity.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert case_data.bond.is_dative() == case_data.is_dative
