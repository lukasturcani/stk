def test_get_vertex2_id(test_case):
    _test_get_vertex2_id(test_case.edge, test_case.vertex2_id)


def _test_get_vertex2_id(edge, id):
    assert edge.get_vertex2_id() == id
