def test_get_id(case_data):
    _test_get_id(case_data.edge, case_data.id)


def _test_get_id(edge, id):
    assert edge.get_id() == id
