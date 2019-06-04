def test_cycle_atoms(cycle_su):
    assert sorted(cycle_su.cycle_atoms()) == list(range(3, 13))
