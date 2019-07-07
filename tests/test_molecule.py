import numpy as np
import stk
import os
from os.path import join


test_dir = 'molecule_tests_output'
if not os.path.exists(test_dir):
    os.mkdir(test_dir)


def test_apply_displacement(tmp_amine2):
    before_conf0 = tmp_amine2.get_position_matrix(conformer_id=0)
    before_conf1 = tmp_amine2.get_position_matrix(conformer_id=1)

    assert np.allclose(before_conf0*4, before_conf1, 1e-4)

    tmp_amine2.apply_displacement([0, 0, 0])

    assert np.allclose(
        a=before_conf0,
        b=tmp_amine2.get_position_matrix(conformer_id=0),
        atol=1e-6
    )
    assert np.allclose(
        a=before_conf1,
        b=tmp_amine2.get_position_matrix(conformer_id=1),
        atol=1e-6
    )

    tmp_amine2.apply_displacement([10, 20, 30])

    after_conf0 = tmp_amine2.get_position_matrix(conformer_id=0)
    assert np.allclose(
        a=before_conf0+[10, 20, 30],
        b=after_conf0,
        atol=1e-6
    )
    assert np.allclose(
        a=before_conf1,
        b=tmp_amine2.get_position_matrix(conformer_id=1),
        atol=1e-6
    )

    tmp_amine2.apply_displacement([-10, 20, -30], conformer_id=1)

    assert np.allclose(
        a=after_conf0,
        b=tmp_amine2.get_position_matrix(conformer_id=0),
        atol=1e-6
    )
    assert np.allclose(
        a=before_conf0+[-10, 20, -30],
        b=tmp_amine2.get_position_matrix(conformer_id=1),
        atol=1e-6
    )


def test_apply_rotation_about_axis():
    assert False


def test_apply_rotation_between_vectors():
    assert False


def apply_rotation_to_minimize_theta():
    assert False


def test_get_atom_coords(tmp_amine2):
    num_atoms = len(tmp_amine2.atoms)
    new_coords = np.stack([np.arange(num_atoms) for _ in range(3)])
    tmp_amine2.set_position_matrix(new_coords)

    for i, atom_coords in enumerate(tmp_amine2.get_atom_coords()):
        assert all(atom_coords == [i, i, i])

    tmp_amine2.set_position_matrix(new_coords*10, conformer_id=1)
    atom_ids = [0, 2, 4]
    coords = tmp_amine2.get_atom_coords(
        atom_ids=atom_ids,
        conformer_id=1
    )
    for atom_id, atom_coords in zip(atom_ids, coords):
        assert all(atom_coords == [atom_id*10]*3)


def test_get_atom_distance(tmp_amine2):
    num_atoms = len(tmp_amine2.atoms)
    coords = np.stack([np.arange(num_atoms) for _ in range(3)])
    tmp_amine2.set_position_matrix(coords)
    for i in range(1, num_atoms):
        assert tmp_amine2.get_atom_distance(i-1, i) == 1
        assert tmp_amine2.get_atom_distance(i, i-1) == 1


def test_get_center_of_mass(tmp_amine2):
    num_atoms = len(tmp_amine2.atoms)
    tmp_amine2.set_position_matrix(np.zeros((num_atoms, 3)))
    assert all(tmp_amine2.get_center_of_mass() == [0, 0, 0])


def test_get_centroid(tmp_amine2):
    num_atoms = len(tmp_amine2.atoms)
    tmp_amine2.set_position_matrix(np.zeros((num_atoms, 3)))
    assert all(tmp_amine2.get_centroid() == [0, 0, 0])

    num_atoms = len(tmp_amine2.atoms)
    tmp_amine2.set_position_matrix(
        position_matrix=np.ones((num_atoms, 3)),
        conformer_id=1
    )
    assert all(tmp_amine2.get_centroid(conformer_id=1) == [1, 1, 1])

    coords = np.stack([np.arange(num_atoms) for _ in range(3)])
    tmp_amine2.set_position_matrix(coords)
    assert np.allclose(
        a=tmp_amine2.get_centroid(atom_ids=[1, 3]),
        b=[2, 2, 2],
        atol=1e-6
    )


def test_get_direction(tmp_amine2):
    direction = tmp_amine2.get_direction(conformer_id=1)
    tmp_amine2.set_orientation(direction, [0, 1, 0], conformer_id=1)
    direction = tmp_amine2.direction(conformer_id=1)
    assert np.allclose(direction, [0, 1, 0], atol=1e-4)


def test_get_maximum_diamter(tmp_amine2):
    # Make a position matrix which sets all atoms to the origin except
    # 1 and 13. These should be placed a distance of 100 apart.
    pos_mat = np.zeros((len(tmp_amine2.atoms), 3))
    pos_mat[1] = [0, -50, 0]
    pos_mat[13] = [0, 50, 0]
    tmp_amine2.set_position_matrix(pos_mat)
    assert abs(tmp_amine2.get_maximum_diameter() - 100) < 1e-6


def test_get_plane_normal():
    assert False


def test_position_matrix(tmp_amine2):
    zeros = np.zeros((len(tmp_amine2.atoms), 3))
    ones = np.ones((len(tmp_amine2.atoms), 3))
    tmp_amine2.set_position_matrix(zeros)
    tmp_amine2.set_position_matrix(ones, conformer_id=1)

    assert np.allclose(zeros, tmp_amine2.get_position_matrix(), 1e-6)
    assert np.allclose(
        a=ones,
        b=tmp_amine2.get_position_matrix(conformer_id=1),
        atol=1e-6
    )


def test_set_centroid(tmp_amine2):
    tmp_amine2.set_centroid([12, 13, 15])
    assert np.allclose(tmp_amine2.get_centroid(), [12, 13, 15], 1e-6)

    tmp_amine2.set_centroid([-12, 4, 160], atom_ids=[1, 3])
    assert not np.allclose(
        a=tmp_amine2.get_centroid(),
        b=[-12, 4, 160],
        atol=1e-6
    )
    assert not np.allclose(
        a=tmp_amine2.get_centroid(atom_ids=[1, 3]),
        b=[-12, 4, 160],
        atol=1e-6
    )


def test_update_from_mae(tmp_amine2, mae_path):
    tmp_amine2.update_from_file(mae_path, None)
    d1 = tmp_amine2.get_maximum_diameter(conformer_id=0)
    d2 = tmp_amine2.get_maximum_diameter(conformer_id=2)
    assert abs(d1 - d2) > 1


def test_update_from_mol(tmp_amine2):
    path = join(test_dir, 'update_from_mol.mol')
    tmp_amine2.write(path=path, conformer_id=1)
    tmp_amine2.update_from_file(path=path, conformer_id=None)

    assert np.allclose(
        a=tmp_amine2.get_position_matrix(conformer_id=1),
        b=tmp_amine2.get_position_matrix(conformer_id=2),
        atol=1e-4
    )
    assert not np.allclose(
        a=tmp_amine2.get_position_matrix(conformer_id=1),
        b=tmp_amine2.get_position_matrix(conformer_id=2),
        atol=1e-4
    )


def test_update_from_xyz(tmp_amine2):
    path = join(test_dir, 'update_from_xyz.xyz')
    tmp_amine2.write(path=path, conformer_id=1)
    tmp_amine2.update_from_file(path=path, conformer_id=None)

    assert np.allclose(
        a=tmp_amine2.get_position_matrix(conformer_id=1),
        b=tmp_amine2.get_position_matrix(conformer_id=2),
        atol=1e-4
    )
    assert not np.allclose(
        a=tmp_amine2.get_position_matrix(conformer_id=1),
        b=tmp_amine2.get_position_matrix(conformer_id=2),
        atol=1e-4
    )


def test_write_pdb(amine2):
    path = join(test_dir, 'test_write.pdb')
    amine2.write(path=path, conformer_id=1)
    bb = stk.BuildingBlock(path)

    assert np.allclose(
        a=amine2.get_position_matrix(conformer_id=1),
        b=bb.get_position_matrix(),
        atol=1e-4
    )
    assert not np.allclose(
        a=amine2.get_position_matrix(conformer_id=1),
        b=bb.get_position_matrix(),
        atol=1e-4
    )
