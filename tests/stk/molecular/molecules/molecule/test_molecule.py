from fixtures import *

import numpy as np
import stk


def test_apply_displacement(molecule, displacement):
    before = molecule.get_position_matrix()
    molecule.apply_displacement(displacement)
    assert np.allclose(
        a=before+displacement,
        b=molecule.get_position_matrix(),
        atol=1e-32,
    )


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

    axis_matrix = np.repeat([axis], len(molecule.atoms), 0).T
    positions = molecule.get_position_matrix() - origin
    return positions - (axis_matrix * (positions @ axis)).T


def test_apply_rotation_about_axis(molecule, angle, axis, origin):
    before = rotational_space_positions(molecule, axis, origin)
    molecule.apply_rotation_about_axis(angle, axis, origin)
    after = rotational_space_positions(molecule, axis, origin)

    for atom_id in range(len(molecule.atoms)):
        applied_rotation = stk.vector_angle(
            vector1=before[atom_id],
            vector2=after[atom_id],
        )
        assert abs(abs(angle) - applied_rotation) < 1e-13


def test_get_atom_positions(get_atom_positions_test_case):
    molecule = get_atom_positions_test_case.molecule
    atom_ids = get_atom_positions_test_case.atom_ids
    true_positions = get_atom_positions_test_case.true_positions

    i = -1
    positions = enumerate(molecule.get_atom_positions(atom_ids))
    for i, position in positions:
        assert np.all(np.equal(position, true_positions[i]))

    assert i+1 == len(atom_ids)
