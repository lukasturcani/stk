import numpy as np
import stk

from ...utilities import is_clone


def test_with_rotation_about_axis(molecule, angle, axis, get_origin):
    """
    Test :meth:`.Molecule.with_rotation_about_axis`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule to test.

    angle : :class:`float`
        The angle of the rotation.

    axis : :class:`numpy.ndarray`
        The axis on which the rotation happens.

    get_origin : :class:`callable`
        Takes a single parameter, `molecule`, and returns a valid
        `origin` parameter for this test. Using a :class:`callable`
        allows the origin for the test to be set to molecule-specific
        values, such as the centroid of the molecule being tested.

    Returns
    -------
    None : :class:`NoneType`

    """

    origin = get_origin(molecule)
    original = rotational_space_positions(molecule, axis, origin)
    new = molecule.with_rotation_about_axis(angle, axis, origin)
    is_clone(new, molecule)
    rotated = rotational_space_positions(new, axis, origin)
    is_rotated(original, rotated, angle)


def is_rotated(original, rotated, angle):
    """
    Check that `rotated` is `original` with a rotation of `angle`.

    Parameters
    ----------
    original : :class:`numpy.ndarray`
        The position matrix of a molecule before rotation.

    rotated : :class:`numpy.ndarray`
        The position matrix of a molecule after rotation.

    angle : :class:`float`
        The rotational angle.

    Returns
    -------
    None : :class:`NoneType`

    """

    for atom_id in range(max(len(original), len(rotated))):
        # No rotation is expected if the atom was located at the
        # origin of the rotation.
        if np.allclose(original[atom_id], [0, 0, 0], 1e-13):
            assert np.allclose(
                a=original[atom_id],
                b=rotated[atom_id],
                atol=1e-13,
            )
        else:
            applied_rotation = stk.vector_angle(
                vector1=original[atom_id],
                vector2=rotated[atom_id],
            )
            assert abs(abs(angle) - applied_rotation) < 1e-13


def rotational_space_positions(molecule, axis, origin):
    """
    Get the atomic coordinates on the plane of the rotation.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule being rotated.

    axis : :class:`numpy.ndarray`
        The axis about which the rotation happens.

    origin : :class:`numpy.ndarray`
        The origin about which the rotation happens.

    Returns
    -------
    :class:`numpy.ndarray`
        An ``[n, 3]`` of atomic positions of `molecule`, projected
        onto the plane about which the rotation happens. The
        `axis` is the normal to this plane.

    """

    axis_matrix = np.repeat([axis], molecule.get_num_atoms(), 0).T
    positions = molecule.get_position_matrix() - origin
    return positions - (axis_matrix * (positions @ axis)).T
