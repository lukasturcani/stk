def test_get_mass(case_data):
    assert case_data.atom.get_mass() == case_data.mass
