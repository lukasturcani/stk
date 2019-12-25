import numpy as np


def test_unchanged(molecule1, molecule2):
    assert np.all(np.equal(
        molecule1.get_position_matrix(),
        molecule2.get_position_matrix(),
    ))
