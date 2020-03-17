def test_get_atom2(case_data):
    assert case_data.bond.get_atom2() is case_data.atom2
