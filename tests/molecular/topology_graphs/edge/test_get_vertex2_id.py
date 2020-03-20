def test_get_vertex2_id(case_data):
    """
    Test :meth:`.Edge.get_vertex2_id`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the edge to test and the correct id of the
        second vertex.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert case_data.edge.get_vertex2_id() == case_data.vertex2_id
