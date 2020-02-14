import numpy as np


def test_with_position_matrix(molecule, get_position_matrix):
    position_matrix = molecule.get_position_matrix()
    _test_with_position_matrix(molecule, get_position_matrix)
    # Test immutability.
    assert np.all(np.equal(
        position_matrix,
        molecule.get_position_matrix(),
    ))


def _test_with_position_matrix(molecule, get_position_matrix):
    position_matrix = get_position_matrix(molecule)
    molecule = molecule.with_position_matrix(position_matrix)
    assert np.all(np.equal(
        position_matrix,
        molecule.get_position_matrix(),
    ))
