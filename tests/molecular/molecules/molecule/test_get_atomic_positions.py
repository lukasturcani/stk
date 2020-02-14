import itertools as it
import numpy as np

from ..utilities import normalize_ids


def test_get_atomic_positions(
    molecule,
    get_position_matrix,
    get_atom_ids,
):
    position_matrix = get_position_matrix(molecule)
    molecule = molecule.with_position_matrix(position_matrix)
    positions = it.zip_longest(
        normalize_ids(molecule, get_atom_ids(molecule)),
        molecule.get_atomic_positions(get_atom_ids(molecule)),
    )
    for atom_id, position in positions:
        assert np.all(np.equal(position, position_matrix[atom_id]))
