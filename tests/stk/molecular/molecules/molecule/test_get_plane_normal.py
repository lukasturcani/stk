import numpy as np
import pytest
import stk

from .utilities import get_num_atom_ids


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
        lambda molecule: range(min(molecule.get_num_atoms(), 3)),
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
        [0., 0., 1],
        [2., -1., 1.32],
    ),
)
def normal(request):
    """
    A plane normal.

    """

    return stk.normalize_vector(np.array(request.param))


def test_get_plane_normal(molecule, get_atom_ids, normal):
    num_atom_ids = get_num_atom_ids(molecule, get_atom_ids)
    if num_atom_ids == 1:
        # Any non-0 vector is valid in this case.
        assert not np.allclose(
            a=[0, 0, 0],
            b=molecule.get_plane_normal(get_atom_ids(molecule)),
            atol=1e-13,
        )
        return

    elif num_atom_ids == 2:
        # Any perpendicular vector is valid in this case.
        result = molecule.get_plane_normal(get_atom_ids(molecule))
        assert result @ get_direction(molecule, 0, 1) < 1e-13
        return

    position_matrix = get_position_matrix(
        molecule=molecule,
        atom_ids=get_atom_ids(molecule),
        normal=normal,
    )
    molecule = molecule.with_position_matrix(position_matrix)
    result = molecule.get_plane_normal(get_atom_ids(molecule))
    # The normal may be parallel or anti-parallel.
    assert abs(abs(result @ normal) - 1) < 1e-13


def get_position_matrix(molecule, atom_ids, normal):
    if atom_ids is None:
        atom_ids = range(molecule.get_num_atoms())
    elif not isinstance(atom_ids, (list, tuple)):
        atom_ids = tuple(atom_ids)

    position_matrix = molecule.get_position_matrix()
    for atom_id in atom_ids:
        remove_component(position_matrix, atom_id, normal)
    return position_matrix


def remove_component(position_matrix, atom_id, normal):
    component_magnitude = position_matrix[atom_id] @ normal
    position_matrix[atom_id, :] -= component_magnitude * normal


def get_direction(molecule, atom1, atom2):
    position_matrix = molecule.get_position_matrix()
    return position_matrix[atom1] - position_matrix[atom2]
