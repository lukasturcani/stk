def test_get_periodicity(case_data):
    """
    Test :meth:`.Bond.get_periodicity`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the bond to test and the correct
        periodicity.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert case_data.bond.get_periodicity() == case_data.periodicity
