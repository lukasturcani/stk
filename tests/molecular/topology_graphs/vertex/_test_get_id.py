def test_get_id(case_data):
    """
    Test :meth:`.Vertex.get_id`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the vertex to test and the correct id.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert case_data.vertex.get_id() == case_data.id
