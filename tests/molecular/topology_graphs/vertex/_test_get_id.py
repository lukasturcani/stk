def test_get_id(case_data):
    _test_get_id(case_data.vertex, case_data.id)


def _test_get_id(vertex, id):
    assert vertex.get_id() == id
