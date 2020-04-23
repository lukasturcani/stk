def test_is_periodic(case_data):
    """
    Test :meth:`.Edge.is_periodic`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the edge to test and the truth about its
        periodicity.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert case_data.edge.is_periodic() == case_data.is_periodic
