def test_cycle_atoms(cycle_su):
    assert set(cycle_su.cycle_atoms()) == {3, 4, 5, 6, 7, 8,
                                           9, 10, 11, 12}


def test_cycle_coords(cycle):
    catoms = cycle.cycle_atoms()
    for i, coords in enumerate(cycle.cycle_coords(), 1):
        assert type(coords[0]) == int
        assert type(coords[1]) == float
        assert type(coords[2]) == float
        assert type(coords[3]) == float
        assert len(coords) == 4
    assert len(catoms) == i
