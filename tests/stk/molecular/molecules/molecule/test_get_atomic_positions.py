import itertools as it
import numpy as np
import pytest


@pytest.fixture(
    params=(
        lambda molecule: np.zeros((molecule.get_num_atoms(), 3)),
        lambda molecule: np.array([
            [i, -i, 10.12*i] for i in range(molecule.get_num_atoms())
        ]),
    ),
)
def get_position_matrix(request):
    """
    A function which returns a valid position matrix for a moleclue.

    The function takes 1 parameter, the :class:`.Molecule` instance
    for which it returns a valid position matrix.

    """

    return request.param


def test_get_atomic_positions(
    molecule,
    get_position_matrix,
    get_atom_ids,
):
    position_matrix = get_position_matrix(molecule)
    molecule = molecule.with_position_matrix(position_matrix)
    atom_ids = get_atom_ids(molecule)
    positions = it.zip_longest(
        get_all_atom_ids(molecule) if atom_ids is None else atom_ids,
        molecule.get_atomic_positions(get_atom_ids(molecule)),
    )
    for atom_id, position in positions:
        assert np.all(np.equal(position, position_matrix[atom_id]))


def get_all_atom_ids(molecule):
    return range(molecule.get_num_atoms())
