def test_get_id(test_case):
    _test_get_id(test_case.vertex, test_case.id)


def _test_get_id(vertex, id):
    assert vertex.get_id() == id
