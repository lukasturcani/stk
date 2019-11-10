import pytest
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


def test_get_atom_positions(molecule, get_atom_ids_no_fail):
    position_matrix = molecule.get_position_matrix()
    atom_ids = get_atom_ids_no_fail(molecule)
    positions = it.zip_longest(
        range(len(molecule.atoms)) if atom_ids is None else atom_ids,
        molecule.get_atom_positions(get_atom_ids_no_fail(molecule)),
    )
    for atom_id, position in positions:
        assert np.all(np.equal(position, position_matrix[atom_id]))


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


def test_get_center_of_mass(molecule, get_atom_ids):
    atom_ids = get_atom_ids(molecule)
    center_of_mass = molecule.get_center_of_mass(atom_ids),

    if atom_ids is None:
        atom_ids = range(len(molecule.atoms))
    else:
        atom_ids = list(get_atom_ids(molecule))

    valid_atom_ids = set(atom_ids)
    atoms = filter(
        lambda atom: atom.id in valid_atom_ids,
        molecule.atoms
    )
    masses = [[atom.mass] for atom in atoms]
    true_center_of_mass = np.divide(
        np.sum(
            a=masses*molecule.get_position_matrix()[atom_ids, :],
            axis=0,
        ),
        sum(mass for mass, in masses),
    )
    assert np.allclose(
        a=true_center_of_mass,
        b=center_of_mass,
        atol=1e-32,
    )


def test_get_centroid(molecule, get_atom_ids):
    atom_ids = get_atom_ids(molecule)
    centroid = molecule.get_centroid(atom_ids)

    if atom_ids is None:
        atom_ids = range(len(molecule.atoms))
    else:
        atom_ids = list(get_atom_ids(molecule))

    true_centroid = np.divide(
        np.sum(
            a=molecule.get_position_matrix()[atom_ids, :],
            axis=0
        ),
        len(atom_ids),
    )
    assert np.allclose(
        a=true_centroid,
        b=centroid,
        atol=1e-32,
    )


class TestGetDirection:

    class _Test1:
        @staticmethod
        def case1():
            bb = stk.BuildingBlock('NCCN')
            bb.set_position_matrix(
                np.array([[i, 0, 0] for i in range(len(bb.atoms))])
            )
            return bb, [1, 0, 0]

    @pytest.mark.parametrize(
        'molecule,direction',
        [
            _Test1.case1()
        ]
    )
    def test_get_direction_1(self, molecule, get_atom_ids, direction):
        atom_ids = get_atom_ids(molecule)
        assert np.allclose(
            a=molecule.get_direction(atom_ids),
            b=direction,
            atol=1e-32,
        )

    class _Test2:
        @staticmethod
        def case1():
            bb = stk.BuildingBlock('NCCN')
            coords = bb.get_position_matrix()
            coords[[1, 3]] = [[1, 1, 1], [3, 3, 3]]
            bb.set_position_matrix(coords)
            return bb, [1/np.sqrt(3)]*3

    @pytest.fixture(params=[
        (1, 3),
    ])
    def atom_ids(self, request):
        return request.param

    @pytest.mark.parametrize(
        'molecule,direction',
        [
            _Test2.case1()
        ]
    )
    def test_get_direction_2(self, molecule, atom_ids, direction):
        assert np.allclose(
            a=molecule.get_direction(atom_ids),
            b=[1/np.sqrt(3), 1/np.sqrt(3), 1/np.sqrt(3)],
            atol=1e-32,
        )
