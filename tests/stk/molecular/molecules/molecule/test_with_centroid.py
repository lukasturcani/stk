import pytest
import numpy as np


@pytest.fixture
def centroid(origin):
    return origin


def test_with_centroid(molecule, get_atom_ids, centroid):
    position_matrix = molecule.get_position_matrix()
    _test_with_centroid(molecule, get_atom_ids, centroid)
    # Test immutability.
    assert np.all(np.equal(
        position_matrix,
        molecule.get_position_matrix(),
    ))


def _test_with_centroid(molecule, get_atom_ids, centroid):
    molecule = molecule.with_centroid(
        position=centroid,
        atom_ids=get_atom_ids(molecule),
    )
    assert np.allclose(
        a=centroid,
        b=molecule.get_centroid(atom_ids=get_atom_ids(molecule)),
        atol=1e-32,
    )
