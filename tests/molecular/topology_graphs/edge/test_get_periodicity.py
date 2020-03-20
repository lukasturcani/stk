def test_get_periodicity(case_data):
    """
    Test :meth:`.Edge.get_periodicity`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the edge to test and the correct
        periodicity.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert case_data.edge.get_periodicity() == case_data.periodicity
