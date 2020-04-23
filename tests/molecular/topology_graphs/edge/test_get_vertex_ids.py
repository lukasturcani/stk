def test_get_vertex_ids(case_data):
    """
    Test :meth:`.Edge.get_vertex_ids`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the edge to test and the correct vertex ids.

    Returns
    -------
    None : :class:`NoneType`

    """

    id1, id2 = case_data.edge.get_vertex_ids()
    assert id1 == case_data.vertex1_id
    assert id2 == case_data.vertex2_id
