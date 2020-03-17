import numpy as np
import pytest

from ..utilities import get_num_atom_ids, normalize_ids


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
        lambda molecule: range(min(molecule.get_num_atoms(), 2)),
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


def test_get_direction(case_data, get_atom_ids):
    _test_get_direction(
        molecule=case_data.molecule,
        position_matrix=case_data.position_matrix,
        get_atom_ids=get_atom_ids,
    )


def _test_get_direction(molecule, position_matrix, get_atom_ids):
    if get_num_atom_ids(molecule, get_atom_ids) == 1:
        # Any non-0 vector is valid in this case.
        assert not np.allclose(
            a=[0, 0, 0],
            b=molecule.get_direction(get_atom_ids(molecule)),
            atol=1e-13,
        )
        return

    direction = get_direction(
        position_matrix=position_matrix,
        atom_ids=normalize_ids(molecule, get_atom_ids(molecule)),
    )
    result = molecule.get_direction(get_atom_ids(molecule))
    # The direction may be parallel or anti-parallel.
    return abs(abs(result @ direction) - 1) < 1e-13


def get_direction(position_matrix, atom_ids):
    atom_positions = position_matrix[atom_ids, :]
    centered_positions = atom_positions - atom_positions.mean(axis=0)
    return np.linalg.svd(centered_positions)[-1][0]


@pytest.mark.parametrize(
    argnames=('position_matrix', 'atom_ids', 'direction'),
    arvalues=(
        (
            np.array([
                [1., 0., 0.],
                [2., 0., 0.],
                [3., 0., 0.],
            ]),
            (0, 1, 2),
            np.array([1., 0., 0.]),
        ),
        (
            np.array([
                [0., 1., 0.],
                [0., 2., 0.],
                [0., 3., 0.],
            ]),
            (0, 1, 2),
            np.array([0., 1., 0.]),
        ),
        (
            np.array([
                [1., 0., 0.],
                [2., 0., 0.],
                [0., 5., 0.],
            ]),
            (0, 1),
            np.array([1., 0., 0.]),
        ),
    ),
)
def test_get_direction_helper(position_matrix, atom_ids, direction):
    """
    Test the :func:`.get_direction` function, defined in this module.

    """

    assert np.allclose(
        a=get_direction(position_matrix, atom_ids),
        b=direction,
        atol=1e-32,
    )
