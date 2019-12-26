import numpy as np
import stk
import pytest


@pytest.fixture(
    params=(
        lambda molecule: None,
        lambda molecule: range(molecule.get_num_atoms()),
        lambda molecule: range(0, molecule.get_num_atoms(), 2),
        lambda molecule: range(0, min(1, molecule.get_num_atoms())),
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


@pytest.fixture(
    params=(
        [1., 0., 0.],
        [0., 1., 0.],
        [0., 1., 0.],
        stk.normalize_vector([1., 1., 1.]),
        stk.normalize_vector([2.231, 0.75, -1.32]),
    )
)
def direction(request):
    """
    A direction vector.

    """

    return np.array(request.param)


def test_get_direction(molecule, get_atom_ids, direction):
    position_matrix = get_position_matrix(
        molecule=molecule,
        atom_ids=get_atom_ids(molecule),
        direction=direction,
    )
    molecule = molecule.with_position_matrix(position_matrix)
    assert np.allclose(
        a=direction,
        b=molecule.get_direction(get_atom_ids(molecule)),
        atol=1e-32,
    )


def get_position_matrix(molecule, atom_ids, direction):
    """
    Create a position matrix with a specific `direction`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule for which the position matrix is created.

    atom_ids : :class:`iterable` of :class:`int`
        The ids of atoms which are to have a specific direction.
        If ``None``, then all atoms are used.

    direction : :class:`numpy.ndarray`
        The desired direction.

    Returns
    -------
    :class:`numpy.ndarray`
        A position matrix for `molecule`.

    """

    if atom_ids is None:
        atom_ids = range(molecule.get_num_atoms())
    elif not isinstance(atom_ids, (tuple, list)):
        atom_ids = tuple(atom_ids)

    position_matrix = molecule.get_position_matrix()
    for atom_id in atom_ids:
        position_matrix[atom_id] = direction*atom_id
    return position_matrix
