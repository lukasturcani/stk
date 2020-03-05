def test_get_vertex1_id(case_data):
    _test_get_vertex1_id(case_data.edge, case_data.vertex1_id)


def _test_get_vertex1_id(edge, id):
    assert edge.get_vertex1_id() == id
