import numpy as np


def test_get_position_matrix(self, molecule, get_position_matrix):
    position_matrix = get_position_matrix(molecule)
    molecule = molecule.with_position_matrix(position_matrix)
    assert np.all(np.equal(
        position_matrix,
        molecule.get_position_matrix(),
    ))
