import numpy as np
import pytest


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

    Retruns
    -------
    :class:`iterable` of :class:`int`
        An `atom_ids` parameter.

    """

    return request.param


def test_get_centroid(case_data, get_atom_ids):
    _test_get_centroid(
        molecule=case_data.molecule,
        position_matrix=case_data.position_matrix,
        get_atom_ids=get_atom_ids,
    )


def _test_get_centroid(molecule, position_matrix, get_atom_ids):
    assert np.allclose(
        a=get_centroid(position_matrix, get_atom_ids(molecule)),
        b=molecule.get_centroid(get_atom_ids(molecule)),
        atol=1e-32,
    )


def get_centroid(position_matrix, atom_ids):
    if atom_ids is None:
        atom_ids = range(len(position_matrix))
    elif isinstance(atom_ids, int):
        atom_ids = (atom_ids, )
    else:
        atom_ids = len(atom_ids)

    return np.sum(position_matrix[atom_ids, :], axis=0) / len(atom_ids)
