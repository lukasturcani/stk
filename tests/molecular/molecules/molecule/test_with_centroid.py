import pytest
import numpy as np


@pytest.fixture(
    params=(
        lambda molecule: None,
        lambda molecule: range(molecule.get_num_atoms()),
        lambda molecule: range(0, molecule.get_num_atoms(), 2),
        lambda molecule: list(
            range(0, min(1, molecule.get_num_atoms()))
        ),
        lambda molecule: tuple(
            range(0, min(1, molecule.get_num_atoms()))
        ),
        lambda molecule: (
            i for i in range(0, min(1, molecule.get_num_atoms()))
        ),
        pytest.param(
            lambda molecule: (),
            marks=pytest.mark.xfail(strict=True, raises=ValueError),
        ),
        lambda molecule: range(min(molecule.get_num_atoms(), 1)),
    ),
)
def get_atom_ids(request):
    """
    Return an atom_ids parameter for a :class:`.Molecule`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule for which `atom_ids` are returned.

    Returns
    -------
    :class:`iterable` of :class:`int`
        An `atom_ids` parameter.

    """

    return request.param


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
        atol=1e-14,
    )
