import itertools as it
import numpy as np

from ..utilities import normalize_ids


def test_get_atomic_positions(case_data, get_atom_ids):
    _test_get_atomic_positions(
        molecule=case_data.molecule,
        position_matrix=case_data.position_matrix,
        get_atom_ids=get_atom_ids,
    )


def _test_get_atomic_positions(
    molecule,
    position_matrix,
    get_atom_ids,
):
    for atom_id, position in it.zip_longest(
        normalize_ids(molecule, get_atom_ids(molecule)),
        molecule.get_atomic_positions(get_atom_ids(molecule)),
    ):
        assert np.all(np.equal(position, position_matrix[atom_id]))
