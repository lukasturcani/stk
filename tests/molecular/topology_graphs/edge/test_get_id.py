def test_get_id(test_case):
    _test_get_id(test_case.edge, test_case.id)


def _test_get_id(edge, id):
    assert edge.get_id() == id
