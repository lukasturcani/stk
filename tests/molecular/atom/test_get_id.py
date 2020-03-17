def test_get_id(case_data):
    assert case_data.atom.get_id() == case_data.id
