def test_get_order(case_data):
    """
    Test :meth:`.Bond.get_order`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the bond to test and the correct bond order.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert case_data.bond.get_order() == case_data.order
