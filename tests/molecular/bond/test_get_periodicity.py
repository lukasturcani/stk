def test_get_periodicity(case_data):
    assert case_data.bond.get_periodicity() == case_data.periodicity
