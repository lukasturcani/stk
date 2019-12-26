import numpy as np
import pytest


def get_random_position_matrix(molecule):
    generator = np.random.RandomState(4)
    return generator.normal(
        loc=23.3,
        scale=32.1,
        size=(molecule.get_num_atoms(), 3),
    )


@pytest.fixture(
    params=(
        lambda molecule: np.zeros((molecule.get_num_atoms(), 3)),
        get_random_position_matrix,
    )
)
def get_position_matrix(request):
    """
    Return a valid position matrix for a molecule.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule for which a position matrix is returned.

    Returns
    -------
    :class:`numpy.ndarray`
        A position matrix for `molecule`.

    """

    return request.param


def test_get_position_matrix(self, molecule, get_position_matrix):
    position_matrix = get_position_matrix(molecule)
    molecule = molecule.with_position_matrix(position_matrix)
    assert np.all(np.equal(
        position_matrix,
        molecule.get_position_matrix(),
    ))
