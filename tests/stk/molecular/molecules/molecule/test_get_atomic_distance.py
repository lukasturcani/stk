import itertools as it
import numpy as np


def test_get_atomic_distance(molecule, get_position_matrix):
    position_matrix = get_position_matrix(molecule)
    molecule = molecule.with_position_matrix(position_matrix)
    check_atomic_distances(molecule, get_distance_matrix(molecule))


def get_distance_matrix(molecule):
    positions_1 = np.repeat(
        a=[molecule.get_position_matrix()],
        repeats=molecule.get_num_atoms(),
        axis=0,
    )
    positions_2 = positions_1.swapaxes(0, 1)
    return np.linalg.norm(positions_1 - positions_2, axis=2)


def check_atomic_distances(molecule, distance_matrix):
    """
    Check that atom distances in `molecule` match `distance_matrix`.

    """

    atom_ids = range(molecule.get_num_atoms())
    for atom1, atom2 in it.product(atom_ids, atom_ids):
        true_distance = distance_matrix[atom1, atom2]
        distance = molecule.get_atomic_distance(atom1, atom2)
        assert abs(true_distance - distance) < 1e-13
