import numpy as np


def test_cycle_atoms(cycle_su):
    assert set(cycle_su.cycle_atoms()) == {3, 4, 5, 6, 7, 8,
                                           9, 10, 11, 12}


def test_cycle_coords(tmp_cycle):
    catoms = tmp_cycle.cycle_atoms()
    new_pos_mat = tmp_cycle.position_matrix(0).T

    # Place cycle atoms at the origin
    for atom in catoms:
        new_pos_mat[atom] = [0 for x in range(3)]

    tmp_cycle.set_position_from_matrix(new_pos_mat.T, 0)

    for i, coords in enumerate(tmp_cycle.cycle_coords(conformer=0), 1):
        assert type(coords[0]) == int
        assert np.allclose(coords[1], np.array([0, 0, 0]), atol=1e-5)
        assert np.allclose(coords[2], np.array([0, 0, 0]), atol=1e-5)
        assert np.allclose(coords[3], np.array([0, 0, 0]), atol=1e-5)
        assert len(coords) == 4
    assert len(catoms) == i
