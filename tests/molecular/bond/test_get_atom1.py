def test_get_atom1(case_data):
    assert case_data.bond.get_atom1() is case_data.atom1
