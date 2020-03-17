def test_get_order(case_data):
    assert case_data.bond.get_order() == case_data.order
