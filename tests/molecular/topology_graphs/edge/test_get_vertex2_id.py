def test_get_vertex2_id(case_data):
    _test_get_vertex2_id(case_data.edge, case_data.vertex2_id)


def _test_get_vertex2_id(edge, id):
    assert edge.get_vertex2_id() == id
