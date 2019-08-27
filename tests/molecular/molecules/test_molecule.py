import numpy as np
import stk
import os
from os.path import join
import pytest


test_dir = 'molecule_tests_output'
if not os.path.exists(test_dir):
    os.mkdir(test_dir)


def test_apply_displacement(tmp_amine2):
    before = tmp_amine2.get_position_matrix()
    tmp_amine2.apply_displacement(np.array([0, 0, 0]))
    assert np.allclose(
        a=before,
        b=tmp_amine2.get_position_matrix(),
        atol=1e-8
    )

    tmp_amine2.apply_displacement(np.array([10, 20, 30]))
    after = tmp_amine2.get_position_matrix()
    assert np.allclose(
        a=before+[10, 20, 30],
        b=after,
        atol=1e-8
    )

    tmp_amine2.apply_displacement(np.array([-10, 20, -30]))
    assert np.allclose(
        a=after+[-10, 20, -30],
        b=tmp_amine2.get_position_matrix(),
        atol=1e-8
    )


def test_apply_rotation_about_axis(tmp_amine2):
    num_atoms = len(tmp_amine2.atoms)
    coords = np.array([[i, 0, 0] for i in range(num_atoms)])
    tmp_amine2.set_position_matrix(coords)

    tmp_amine2.apply_rotation_about_axis(
        angle=-np.pi/2,
        axis=np.array([0, 1, 0]),
        origin=np.array([0, 0, 0])
    )

    assert np.allclose(
        a=tmp_amine2.get_position_matrix(),
        b=np.array([[0, 0, i] for i in range(num_atoms)]),
        atol=1e-8
    )


def test_apply_rotation_between_vectors(tmp_amine2):
    assert not np.allclose(
        a=next(tmp_amine2.get_bonder_direction_vectors())[-1],
        b=[1, 0, 0],
        atol=1e-8
    )

    tmp_amine2.apply_rotation_between_vectors(
        start=next(tmp_amine2.get_bonder_direction_vectors())[-1],
        target=[1, 0, 0],
        origin=tmp_amine2.get_centroid()
    )
    assert np.allclose(
        a=stk.normalize_vector(
            next(tmp_amine2.get_bonder_direction_vectors())[-1]
        ),
        b=[1, 0, 0],
        atol=1e-8
    )


def test_apply_rotation_to_minimize_angle(tmp_amine2):
    num_atoms = len(tmp_amine2.atoms)
    coords = np.array([[i, 0, 0] for i in range(num_atoms)])
    tmp_amine2.set_position_matrix(coords)

    tmp_amine2.apply_rotation_to_minimize_angle(
        start=np.array([1, 0, 0]),
        target=np.array([0, 0, 1]),
        axis=np.array([0, 1, 0]),
        origin=np.array([0, 0, 0])
    )

    assert np.allclose(
        a=tmp_amine2.get_position_matrix(),
        b=np.array([[0, 0, i] for i in range(num_atoms)]),
        atol=1e-8
    )


def test_get_atom_coords(tmp_amine2):
    num_atoms = len(tmp_amine2.atoms)
    new_coords = np.array([[i]*3 for i in range(num_atoms)])
    tmp_amine2.set_position_matrix(new_coords)

    for i, atom_coords in enumerate(tmp_amine2.get_atom_coords()):
        assert all(atom_coords == [i, i, i])

    tmp_amine2.set_position_matrix(new_coords*10)

    # Test different input types.
    all_atom_ids = ((0, 2, 4), [0, 2, 4], (i for i in [0, 2, 4]))
    for atom_ids in all_atom_ids:
        coords = zip(
            (0, 2, 4),
            tmp_amine2.get_atom_coords(atom_ids=atom_ids)
        )
        for i, (atom_id, atom_coords) in enumerate(coords, 1):
            assert all(atom_coords == [atom_id*10]*3)
        assert i == 3


def test_get_atom_distance(tmp_amine2):
    num_atoms = len(tmp_amine2.atoms)
    coords = np.array([[i, 0, 0] for i in range(num_atoms)])
    tmp_amine2.set_position_matrix(coords)
    for i in range(1, num_atoms):
        assert tmp_amine2.get_atom_distance(i-1, i) == 1
        assert tmp_amine2.get_atom_distance(i, i-1) == 1


def test_get_cached_mol(tmp_amine2, aldehyde2):
    try:
        stk.BuildingBlock._cache = {}
        tmp_amine2.update_cache()
        cached = stk.BuildingBlock.get_cached_mol(
            identity_key=tmp_amine2.get_identity_key()
        )
        assert cached is tmp_amine2

        with pytest.raises(KeyError):
            stk.BuildingBlock.get_cached_mol(
                identity_key=aldehyde2.get_identity_key()
            )

        default = object()
        cached = stk.BuildingBlock.get_cached_mol(
            identity_key=aldehyde2.get_identity_key(),
            default=default
        )
        assert default is cached

    except Exception:
        raise

    finally:
        stk.BuildingBlock._cache = {}


def test_get_center_of_mass(tmp_amine2):
    num_atoms = len(tmp_amine2.atoms)
    tmp_amine2.set_position_matrix(np.zeros((num_atoms, 3)))
    assert all(tmp_amine2.get_center_of_mass() == [0, 0, 0])

    new_coords = tmp_amine2.get_position_matrix()
    new_coords[(0, 2, 4), :] = np.ones((3, 3))
    tmp_amine2.set_position_matrix(new_coords)
    assert not all(tmp_amine2.get_center_of_mass() == [1, 1, 1])

    all_atom_ids = ((0, 2, 4), [0, 2, 4], (i for i in [0, 2, 4]))
    for atom_ids in all_atom_ids:
        assert all(
            tmp_amine2.get_center_of_mass(atom_ids) == [1, 1, 1]
        )


def test_get_centroid(tmp_amine2):
    num_atoms = len(tmp_amine2.atoms)
    tmp_amine2.set_position_matrix(np.zeros((num_atoms, 3)))
    assert all(tmp_amine2.get_centroid() == [0, 0, 0])

    num_atoms = len(tmp_amine2.atoms)
    tmp_amine2.set_position_matrix(np.ones((num_atoms, 3)))
    assert all(tmp_amine2.get_centroid() == [1, 1, 1])

    coords = np.array([[i]*3 for i in range(num_atoms)])
    tmp_amine2.set_position_matrix(coords)

    all_atom_ids = ((1, 3), [1, 3], (i for i in [1, 3]))
    for atom_ids in all_atom_ids:
        assert np.allclose(
            a=tmp_amine2.get_centroid(atom_ids=atom_ids),
            b=[2, 2, 2],
            atol=1e-8
        )


def test_get_direction(tmp_amine2):
    num_atoms = len(tmp_amine2.atoms)
    coords = np.array([[i, 0, 0] for i in range(num_atoms)])
    tmp_amine2.set_position_matrix(coords)

    assert np.allclose(
        a=tmp_amine2.get_direction(),
        b=[1, 0, 0],
        atol=1e-8
    )

    coords[[1, 3]] = [[1, 1, 1], [3, 3, 3]]
    tmp_amine2.set_position_matrix(coords)

    all_atom_ids = ((1, 3), [1, 3], (i for i in [1, 3]))
    for atom_ids in all_atom_ids:
        assert np.allclose(
            a=tmp_amine2.get_direction(atom_ids=atom_ids),
            b=stk.normalize_vector([1, 1, 1]),
            atol=1e-8
        )


def test_get_maximum_diamter(tmp_amine2):
    # Make a position matrix which sets all atoms to the origin except
    # 1 and 13. These should be placed a distance of 100 apart.
    pos_mat = np.zeros((len(tmp_amine2.atoms), 3))
    pos_mat[1] = [0, -50, 0]
    pos_mat[13] = [0, 50, 0]
    tmp_amine2.set_position_matrix(pos_mat)
    assert abs(tmp_amine2.get_maximum_diameter() - 100) < 1e-8

    all_atom_ids = (
        [i for i in range(len(tmp_amine2.atoms)) if i not in {1, 13}],
        (i for i in range(len(tmp_amine2.atoms)) if i not in {1, 13}),
        tuple(
            i for i in range(len(tmp_amine2.atoms)) if i not in {1, 13}
        )
    )
    for atom_ids in all_atom_ids:
        assert abs(tmp_amine2.get_maximum_diameter(atom_ids)) < 1e-8


def test_get_plane_normal(tmp_amine2):
    coords = tmp_amine2.get_position_matrix()
    all_atom_ids = ([1, 13], (1, 13), (i for i in [1, 13]))
    coords[[1, 13], 2] = 0
    tmp_amine2.set_position_matrix(coords)

    assert not np.allclose(
        a=tmp_amine2.get_plane_normal(),
        b=[0, 0, 1],
        atol=1e-8
    )

    for atom_ids in all_atom_ids:
        assert np.allclose(
            a=tmp_amine2.get_plane_normal(atom_ids),
            b=[0, 0, 1],
            atol=1e-8
        )

    coords[:, 2] = 0
    tmp_amine2.set_position_matrix(coords)
    assert np.allclose(
        a=tmp_amine2.get_plane_normal(),
        b=[0, 0, 1],
        atol=1e-8
    )


def test_get_set_position_matrix(tmp_amine2):
    zeros = np.zeros((len(tmp_amine2.atoms), 3))
    tmp_amine2.set_position_matrix(zeros)
    assert np.allclose(zeros, tmp_amine2.get_position_matrix(), 1e-8)

    ones = np.ones((len(tmp_amine2.atoms), 3))
    tmp_amine2.set_position_matrix(ones)
    assert np.allclose(ones, tmp_amine2.get_position_matrix(), 1e-8)


def test_has_cached_mol(tmp_amine2, aldehyde2):
    try:
        stk.BuildingBlock._cache = {}
        tmp_amine2.update_cache()
        is_cached = stk.BuildingBlock.has_cached_mol(
            identity_key=tmp_amine2.get_identity_key()
        )
        assert is_cached

        is_cached = stk.BuildingBlock.has_cached_mol(
            identity_key=aldehyde2.get_identity_key()
        )
        assert not is_cached

    except Exception:
        raise

    finally:
        stk.BuildingBlock._cache = {}


def test_set_centroid(tmp_amine2):
    tmp_amine2.set_centroid([12, 13, 15])
    assert np.allclose(tmp_amine2.get_centroid(), [12, 13, 15], 1e-8)

    set_all_atom_ids = ([1, 3], (2, 3), (i for i in [0, 3]))
    get_all_atom_ids = ([1, 3], (2, 3), (i for i in [0, 3]))
    all_atom_ids = zip(set_all_atom_ids, get_all_atom_ids)
    for set_atom_ids, get_atom_ids in all_atom_ids:
        tmp_amine2.set_centroid([-12, 4, 160], atom_ids=set_atom_ids)
        assert not np.allclose(
            a=tmp_amine2.get_centroid(),
            b=[-12, 4, 160],
            atol=1e-8
        )
        assert np.allclose(
            a=tmp_amine2.get_centroid(atom_ids=get_atom_ids),
            b=[-12, 4, 160],
            atol=1e-8
        )


def test_update_from_rdkit_mol(tmp_amine2):
    before = tmp_amine2.get_position_matrix()

    mol = tmp_amine2.to_rdkit_mol()
    conf = mol.GetConformer()
    for atom_id, coord in enumerate(conf.GetPositions()):
        conf.SetAtomPosition(atom_id, 0.5*coord)

    tmp_amine2.update_from_rdkit_mol(mol)
    after = tmp_amine2.get_position_matrix()
    assert np.allclose(conf.GetPositions(), after, 1e-8)
    assert not np.allclose(before, after, 1e-8)


def test_update_from_mae(tmp_amine2, mae_path):
    before = tmp_amine2.get_maximum_diameter()
    tmp_amine2.update_from_file(mae_path)
    after = tmp_amine2.get_maximum_diameter()
    assert abs(before - after) > 1


def test_update_from_mol(tmp_amine2, amine2_conf1):
    assert not np.allclose(
        a=tmp_amine2.get_position_matrix(),
        b=amine2_conf1.get_position_matrix(),
        atol=1e-4
    )

    path = join(test_dir, 'update_from_mol.mol')
    amine2_conf1.write(path=path)
    tmp_amine2.update_from_file(path=path)

    assert np.allclose(
        a=tmp_amine2.get_position_matrix(),
        b=amine2_conf1.get_position_matrix(),
        atol=1e-4
    )


def test_update_from_xyz(tmp_amine2, amine2_conf1):
    assert not np.allclose(
        a=tmp_amine2.get_position_matrix(),
        b=amine2_conf1.get_position_matrix(),
        atol=1e-4
    )

    path = join(test_dir, 'update_from_mol.xyz')
    amine2_conf1.write(path=path)
    tmp_amine2.update_from_file(path=path)

    assert np.allclose(
        a=tmp_amine2.get_position_matrix(),
        b=amine2_conf1.get_position_matrix(),
        atol=1e-4
    )


def test_write_pdb(amine2):
    path = join(test_dir, 'test_write.pdb')
    amine2.write(path=path)
    bb = stk.BuildingBlock.init_from_file(path)

    assert np.allclose(
        a=amine2.get_position_matrix(),
        b=bb.get_position_matrix(),
        atol=1e-4
    )


def test_to_rdkit_mol(tmp_amine2):
    tmp_amine2.atoms[0].charge = 12
    mol = tmp_amine2.to_rdkit_mol()
    assert mol.GetNumConformers() == 1
    assert mol.GetAtomWithIdx(0).GetFormalCharge() == 12
    assert mol.GetNumAtoms() == len(tmp_amine2.atoms)
    for atom, rdkit_atom in zip(tmp_amine2.atoms, mol.GetAtoms()):
        assert atom.atomic_number == rdkit_atom.GetAtomicNum()
        assert atom.mass == rdkit_atom.GetMass()

    assert mol.GetNumBonds() == len(tmp_amine2.bonds)
    for bond, rdkit_bond in zip(tmp_amine2.bonds, mol.GetBonds()):
        assert rdkit_bond.GetBondTypeAsDouble() == bond.order
        assert rdkit_bond.GetBeginAtomIdx() == bond.atom1.id
        assert rdkit_bond.GetEndAtomIdx() == bond.atom2.id


def test_update_cache(tmp_amine2):
    try:
        cache = dict(stk.BuildingBlock._cache)
        # Create a cached molecule.
        cached = stk.BuildingBlock.init_from_rdkit_mol(
            mol=tmp_amine2.to_rdkit_mol(),
            functional_groups=['amine'],
            use_cache=True
        )
        assert tmp_amine2 is not cached
        assert (
            tmp_amine2.get_identity_key() == cached.get_identity_key()
        )
        assert not hasattr(cached, 'test_value')

        tmp_amine2.test_value = 15
        tmp_amine2.update_cache()
        assert cached.test_value == tmp_amine2.test_value

    except Exception:
        raise

    finally:
        stk.BuildingBlock._cache = cache
