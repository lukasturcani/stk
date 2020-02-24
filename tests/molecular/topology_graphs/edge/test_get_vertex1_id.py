def test_get_vertex1_id(test_case):
    _test_get_vertex1_id(test_case.edge, test_case.vertex1_id)


def _test_get_vertex1_id(edge, id):
    assert edge.get_vertex1_id() == id
