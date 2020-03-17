def test_get_charge(case_data):
    assert case_data.atom.get_charge() == case_data.charge
