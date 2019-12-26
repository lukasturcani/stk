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


@pytest.fixture
def centroid(origin):
    return origin


def test_get_centroid(molecule, get_atom_ids, centroid):
    position_matrix = get_position_matrix(
        molecule=molecule,
        atom_ids=get_atom_ids(molecule),
        centroid=centroid,
    )
    molecule = molecule.with_position_matrix(position_matrix)
    assert np.allclose(
        a=centroid,
        b=molecule.get_centroid(get_atom_ids(molecule)),
        atol=1e-32,
    )


def get_position_matrix(molecule, atom_ids, centroid):
    """
    Create a position matrix with a specific `centroid`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule for which the position matrix is created.

    atom_ids : :class:`iterable` of :class:`int`
        The ids of atoms which are to have a specific `centroid`.
        If ``None``, then all atoms are used.

    centroid : :class:`numpy.ndarray`
        The desired centroid.

    Returns
    -------
    :class:`numpy.ndarray`
        A position matrix for `molecule`.


    """

    if atom_ids is None:
        atom_ids = range(molecule.get_num_atoms())
    elif not isinstance(atom_ids, (list, tuple)):
        atom_ids = tuple(atom_ids)

    position_matrix = molecule.get_position_matrix()
    # First, place all atoms on the centroid.
    position_matrix[atom_ids, :] = centroid
    # Displace pairs of atoms such that the center of mass remains
    # same.
    generator = np.random.RandomState(4)
    for atom_id in range(0, molecule.get_num_atoms()-1, 2):
        displacement = generator.normal(scale=59, size=3)
        position_matrix[atom_id] += displacement
        position_matrix[atom_id] -= displacement
    return position_matrix
