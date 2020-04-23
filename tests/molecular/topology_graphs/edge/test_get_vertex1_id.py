def test_get_vertex1_id(case_data):
    """
    Test :meth:`.Edge.get_vertex1_id`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the edge to test and the correct id of the
        first vertex.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert case_data.edge.get_vertex1_id() == case_data.vertex1_id
