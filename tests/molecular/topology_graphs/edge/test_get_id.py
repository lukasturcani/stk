def test_get_id(case_data):
    """
    Test :meth:`.Edge.get_id`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the edge to test and the correct id.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert case_data.edge.get_id() == case_data.id
