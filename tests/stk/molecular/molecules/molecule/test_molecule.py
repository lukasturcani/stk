from fixtures import *

import numpy as np
import itertools as it
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


def test_get_atom_positions(molecule, get_atom_ids):
    position_matrix = molecule.get_position_matrix()
    atom_ids = get_atom_ids(molecule)

    i = -1
    positions = enumerate(molecule.get_atom_positions(atom_ids))
    for i, position in positions:
        atom_id = atom_ids[i]
        assert np.all(np.equal(position, position_matrix[atom_id]))

    assert i+1 == len(atom_ids)


def test_get_atom_distance(molecule):
    position_matrix = molecule.get_position_matrix()
    positions_1 = np.repeat([position_matrix], len(molecule.atoms), 0)
    positions_2 = positions_1.swapaxes(0, 1)
    distance_matrix = np.linalg.norm(positions_1 - positions_2, axis=2)

    atom_ids = range(len(molecule.atoms))
    for atom1, atom2 in it.product(atom_ids, atom_ids):
        true_distance = distance_matrix[atom1, atom2]
        distance = molecule.get_atom_distance(atom1, atom2)
        assert abs(true_distance - distance) < 1e-13


def test_get_centroid(molecule, get_atom_ids):
    atom_ids = get_atom_ids(molecule)
    if len(atom_ids) == 0:
        assert np.all(np.isnan(molecule.get_centroid(atom_ids)))
    else:
        true_centroid = np.divide(
            np.sum(
                a=molecule.get_position_matrix()[atom_ids, :],
                axis=0
            ),
            len(atom_ids),
        )
        centroid = molecule.get_centroid(atom_ids)
        assert np.allclose(
            a=true_centroid,
            b=centroid,
            atol=1e-32,
        )
