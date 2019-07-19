import os
import numpy as np
import stk
import itertools as it

if not os.path.exists('building_block_tests_output'):
    os.mkdir('building_block_tests_output')


def test_init():
    assert False


def test_init_from_random_file():
    assert False


def test_init_from_smiles():
    mol0 = stk.BuildingBlock.init_from_smiles('NCCCN', ['amine'])
    mol1 = stk.BuildingBlock.init_from_smiles('NCCCN', ['amine'])
    assert mol0 is not mol1

    mol2 = stk.BuildingBlock.init_from_smiles(
        smiles='NCCCN',
        functional_groups=['amine'],
        use_cache=True
    )
    mol3 = stk.BuildingBlock.init_from_smiles(
        smiles='NCCCN',
        functional_groups=['amine'],
        use_cache=True
    )
    assert mol0 is not mol2 and mol1 is not mol2
    assert mol2 is mol3

    mol4 = stk.BuildingBlock.init_from_smiles(
        smiles='NCCCN',
        functional_groups=['aldehyde'],
        use_cache=True
    )
    assert mol3 is not mol4


def test_get_bonder_ids(aldehyde3):
    # Make sure that by default all bonder ids are yielded.
    all_ids = []
    for func_group in aldehyde3.func_groups:
        all_ids.extend(func_group.bonder_ids)
    all_ids.sort()
    default_ids = sorted(aldehyde3.get_bonder_ids())
    assert default_ids == all_ids

    # Make sure that providing all fg ids explicitly is the same
    # as default behaviour.
    fg_ids = range(len(aldehyde3.func_groups))
    explicit_ids = sorted(aldehyde3.get_bonder_ids(fg_ids=fg_ids))
    assert default_ids == explicit_ids

    # Make sure when providing a subset of fg ids, only those are
    # returned.
    subset_ids = []
    fgs = [aldehyde3.func_groups[0], aldehyde3.func_groups[2]]
    for func_group in fgs:
        subset_ids.extend(func_group.bonder_ids)
    subset_ids.sort()

    returned_subset_ids = sorted(
        aldehyde3.get_bonder_ids(fg_ids=[0, 2])
    )
    assert returned_subset_ids == subset_ids


def test_get_bonder_centroids(tmp_aldehyde3):
    # Set the position of all bonder atoms to (0, 0, 0).
    bonder_ids = list(tmp_aldehyde3.get_bonder_ids())
    coords = tmp_aldehyde3.get_position_matrix()
    coords[bonder_ids, :] = np.zeros((len(bonder_ids), 3))
    tmp_aldehyde3.set_position_matrix(coords)
    # Check that the bonder centroids are all at (0, 0, 0).
    for i, centroid in enumerate(tmp_aldehyde3.get_bonder_centroids()):
        assert np.allclose(centroid, [0, 0, 0], 1e-6)
    assert i == 2

    # Set the position of the bonder atoms in functional groups 1 and 2
    # to (1, 1, 1).
    fg_ids = [1, 2]
    bonder_ids = list(tmp_aldehyde3.get_bonder_ids(fg_ids=fg_ids))
    coords[bonder_ids, :] = np.ones((len(bonder_ids), 3))
    tmp_aldehyde3.set_position_matrix(coords)
    # Check that the bonder centroids of functional groups 1 and 2 are
    # at (1, 1, 1).
    centroids = tmp_aldehyde3.get_bonder_centroids(fg_ids=[1, 2])
    for i, centroid in enumerate(centroids):
        assert np.allclose(centroid, [1, 1, 1], 1e-6)
    assert i == 1
    # Check that the bonder centroid of functional group 0 is still at
    # (0, 0, 0).
    for i, centroid in tmp_aldehyde3.get_bonder_centroids(fg_ids=[0]):
        assert np.allclose(centroid, [0, 0, 0], 1e-6)
    assert i == 0


def test_get_bonder_plane(tmp_amine4):
    # First check that when 3 fgs are used, the bonder centroids all
    # sit on the plane.
    for fg_ids in it.combinations(range(4), 3):
        a, b, c, d = tmp_amine4.get_bonder_plane(fg_ids=fg_ids)
        for x, y, z in tmp_amine4.get_bonder_centroids(fg_ids=fg_ids):
            product = a*x + b*y + c*z
            assert abs(product-d) < 1e-6

    # When 4 are used make sure that a best fit plane is produced.
    # Ensure that centroids are placed such that best fit plane is a
    # distance of 0.5 away from all centroids.
    coords = tmp_amine4.get_position_matrix()
    bonder_ids = list(tmp_amine4.get_bonder_ids())
    coords[bonder_ids[0]] = [np.sqrt(2), 0, 0]
    coords[bonder_ids[1]] = [0, 0, np.sqrt(2)]
    coords[bonder_ids[2]] = [1, -1, 0]
    coords[bonder_ids[3]] = [1, 1, 0]
    tmp_amine4.set_position_matrix(coords)

    a, b, c, d = tmp_amine4.get_bonder_plane()
    for x, y, z in tmp_amine4.get_bonder_centroids():
        product = a*x + b*y + c*z
        assert np.allclose(abs(product-d), 0.5, 1e-6)


def test_get_bonder_plane_normal(tmp_amine2):
    coords = tmp_amine2.get_position_matrix()
    coords[:, 2] = 0
    tmp_amine2.set_position_matrix(coords)
    assert np.allclose(
        a=tmp_amine2.get_plane_normal(),
        b=[0, 0, 1],
        atol=1e-6
    )


def test_get_bonder_distances(tmp_amine2):
    coords = tmp_amine2.get_position_matrix()
    atoms0 = tmp_amine2.func_groups[0].bonder_ids
    coords[atoms0, :] = np.zeros((len(atoms0), 3))

    atoms1 = tmp_amine2.func_groups[1].bonder_ids
    coords[atoms1, :] = np.zeros((len(atoms1), 3))
    coords[atoms1, 0] = np.ones((len(atoms1, )))

    for fg1, fg2, distance in tmp_amine2.get_bonder_distances():
        assert abs(distance - 1) < 1e-6


def test_get_bonder_direction_vectors(tmp_aldehyde3):
    pos_mat = tmp_aldehyde3.get_position_matrix()
    # Set the coordinate of each bonder to the id of the fg.
    for fg in tmp_aldehyde3.func_groups:
        for bonder in fg.bonder_ids:
            pos_mat[bonder, :] = [fg.id, fg.id, fg.id]
    tmp_aldehyde3.set_position_matrix(pos_mat)

    dir_vectors = tmp_aldehyde3.get_bonder_direction_vectors()
    for i, (id1, id2, v) in enumerate(dir_vectors):
        # Calculate the expected direction vector based on ids.
        d = stk.normalize_vector(np.array([id2]*3) - np.array([id1]*3))
        assert np.allclose(d, v, atol=1e-8)
    assert i == 2


def test_get_centroid_centroid_direction_vector(aldehyde3):
    c1 = aldehyde3.get_centroid(atom_ids=aldehyde3.get_bonder_ids())
    c2 = aldehyde3.get_centroid()
    assert np.allclose(
        a=stk.normalize_vector(c2-c1),
        b=aldehyde3.get_centroid_centroid_direction_vector(),
        atol=1e-8
    )


def test_get_functional_groups(amine2):
    amines = amine2.get_functional_groups(['amine'])

    assert amines[0].atom_ids == (0, 5, 6)
    assert amines[0].bonder_ids == (0, )
    assert amines[0].deleter_ids == (5, 6)

    assert amines[1].atom_ids == (4, 13, 14)
    assert amines[1].bonder_ids == (4, )
    assert amines[1].deleter_ids == (13, 14)

    aldehydes = amine2.get_functional_groups(['aldehyde'])
    assert not aldehydes


def test_dump_and_load(tmp_amine2):
    path = os.path.join('building_block_tests_output', 'mol.dump')

    tmp_amine2.test_attr1 = 'something'
    tmp_amine2.test_attr2 = 12
    tmp_amine2.test_attr3 = ['12', 'something', 21]
    tmp_amine2.test_attr4 = 'skip'
    include_attrs = ['test_attr1', 'test_attr2', 'test_attr3']

    # Add some custom atom properties.
    tmp_amine2.atoms[0].some_prop = 'custom atom prop'

    tmp_amine2.dump(path, include_attrs)
    mol2 = stk.Molecule.load(path)

    assert tmp_amine2 is not mol2
    assert mol2.func_groups == tmp_amine2.func_groups

    assert tmp_amine2.test_attr1 == mol2.test_attr1
    assert tmp_amine2.test_attr2 == mol2.test_attr2
    assert tmp_amine2.test_attr3 == mol2.test_attr3
    assert not hasattr(mol2, 'test_attr4')
    assert vars(tmp_amine2.atoms[0]) == vars(mol2.atoms[0])
    assert vars(tmp_amine2.atoms[1]) == vars(mol2.atoms[1])

    mol3 = stk.Molecule.load(path, use_cache=True)
    assert mol3 is not mol2
    mol4 = stk.Molecule.load(path, use_cache=True)
    assert mol3 is mol4


def test_shift_fgs(amine4):
    ids = [10, 20, 30, 40]
    shifted = amine4.shift_fgs(ids, 32)

    for i, (fg1, fg2) in enumerate(zip(amine4.func_groups, shifted)):
        assert fg1 is not fg2
        assert fg2.id == ids[i]

        for a1, a2 in zip(fg1.atom_ids, fg2.atom_ids):
            assert a1 + 32 == a2

        for a1, a2 in zip(fg1.bonder_ids, fg2.bonder_ids):
            assert a1 + 32 == a2

        for a1, a2 in zip(fg1.deleter_ids, fg2.deleter_ids):
            assert a1 + 32 == a2
