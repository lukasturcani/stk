import os
import numpy as np
import stk
import itertools as it
from collections import Counter
from os.path import join
import rdkit.Chem.AllChem as rdkit


if not os.path.exists('building_block_tests_output'):
    os.mkdir('building_block_tests_output')


def test_init_rdkit():
    rdkit_mol = rdkit.AddHs(rdkit.MolFromSmiles('NCCCN'))
    rdkit.EmbedMolecule(rdkit_mol, rdkit.ETKDGv2())

    mol0 = stk.BuildingBlock(rdkit_mol, ['amine'])
    # Test that all values are initialized correctly.
    assert len(mol0.func_groups) == 2
    fg_types = stk.dedupe(fg.info.name for fg in mol0.func_groups)
    assert sum(1 for _ in fg_types) == 1
    assert len(mol0.atoms) == 15
    assert len(mol0.bonds) == 14

    atom_count = {
        (stk.H, 0): 10,
        (stk.N, 0): 2,
        (stk.C, 0): 3
    }
    assert atom_count == Counter(
        (a.__class__, a.charge) for a in mol0.atoms
    )

    expected_bonds = {
        frozenset({stk.N, stk.C}): 2,
        frozenset({stk.C}): 2,
        frozenset({stk.H, stk.N}): 4,
        frozenset({stk.H, stk.C}): 6
    }
    assert expected_bonds == Counter(
        frozenset({b.atom1.__class__, b.atom2.__class__})
        for b in mol0.bonds
    )

    # Test that caching is working properly.
    mol1 = stk.BuildingBlock(rdkit_mol, ['amine'])
    assert mol0 is not mol1

    mol2 = stk.BuildingBlock(
        mol=rdkit_mol,
        functional_groups=['amine'],
        use_cache=True
    )
    mol3 = stk.BuildingBlock(
        mol=rdkit_mol,
        functional_groups=['amine'],
        use_cache=True
    )
    assert mol0 is not mol2 and mol1 is not mol2
    assert mol2 is mol3

    mol4 = stk.BuildingBlock(
        mol=rdkit_mol,
        functional_groups=['aldehyde'],
        use_cache=True
    )
    assert mol3 is not mol4

    # Make sure that charged molecules are handled correctly.
    negative_carbon = rdkit.AddHs(rdkit.MolFromSmiles('NC[C-]CN'))
    rdkit.EmbedMolecule(negative_carbon, rdkit.ETKDGv2())

    mol5 = stk.BuildingBlock(
        mol=negative_carbon,
        functional_groups=['amine'],
        use_cache=True
    )
    assert mol5 is not mol0
    # Test that all values are initialized correctly.
    assert len(mol5.func_groups) == 2
    fg_types = stk.dedupe(fg.info.name for fg in mol5.func_groups)
    assert sum(1 for _ in fg_types) == 1
    assert len(mol5.atoms) == 13
    assert len(mol5.bonds) == 12

    atom_count = {
        (stk.C, 0): 2,
        (stk.C, -1): 1,
        (stk.N, 0): 2,
        (stk.H, 0): 8,
    }
    assert atom_count == Counter(
        (a.__class__, a.charge) for a in mol5.atoms
    )

    expected_bonds = {
        frozenset({stk.N, stk.C}): 2,
        frozenset({stk.C}): 2,
        frozenset({stk.H, stk.N}): 4,
        frozenset({stk.H, stk.C}): 4
    }
    assert expected_bonds == Counter(
        frozenset({b.atom1.__class__, b.atom2.__class__})
        for b in mol5.bonds
    )

    negative_nitrogen = rdkit.AddHs(rdkit.MolFromSmiles('[N-]CCCN'))
    rdkit.EmbedMolecule(negative_nitrogen, rdkit.ETKDGv2())

    mol6 = stk.BuildingBlock(
        mol=negative_nitrogen,
        functional_groups=['amine'],
        use_cache=True
    )
    assert mol6 is not mol5 and mol6 is not mol0
    # Test that all values are initialized correctly.
    assert len(mol6.func_groups) == 1
    fg_types = stk.dedupe(fg.info.name for fg in mol6.func_groups)
    assert sum(1 for _ in fg_types) == 1
    assert len(mol6.atoms) == 13
    assert len(mol6.bonds) == 12

    atom_count = {
        (stk.C, 0): 3,
        (stk.N, 0): 1,
        (stk.N, -1): 1,
        (stk.H, 0): 8,
    }
    assert atom_count == Counter(
        (a.__class__, a.charge) for a in mol6.atoms
    )

    expected_bonds = {
        frozenset({stk.N, stk.C}): 2,
        frozenset({stk.C}): 2,
        frozenset({stk.H, stk.N}): 2,
        frozenset({stk.H, stk.C}): 6
    }
    assert expected_bonds == Counter(
        frozenset({b.atom1.__class__, b.atom2.__class__})
        for b in mol6.bonds
    )


def test_init_mol(bb_dir):
    mol0 = stk.BuildingBlock(join(bb_dir, 'neutral.mol'), ['amine'])
    # Test that all values are initialized correctly.
    assert len(mol0.func_groups) == 2
    fg_types = stk.dedupe(fg.info.name for fg in mol0.func_groups)
    assert sum(1 for _ in fg_types) == 1
    assert len(mol0.atoms) == 15
    assert len(mol0.bonds) == 14

    atom_count = {
        (stk.H, 0): 10,
        (stk.N, 0): 2,
        (stk.C, 0): 3
    }
    assert atom_count == Counter(
        (a.__class__, a.charge) for a in mol0.atoms
    )

    expected_bonds = {
        frozenset({stk.N, stk.C}): 2,
        frozenset({stk.C}): 2,
        frozenset({stk.H, stk.N}): 4,
        frozenset({stk.H, stk.C}): 6
    }
    assert expected_bonds == Counter(
        frozenset({b.atom1.__class__, b.atom2.__class__})
        for b in mol0.bonds
    )

    # Test that caching is working properly.
    mol1 = stk.BuildingBlock(join(bb_dir, 'neutral.mol'), ['amine'])
    assert mol0 is not mol1

    mol2 = stk.BuildingBlock(
        mol=join(bb_dir, 'neutral.mol'),
        functional_groups=['amine'],
        use_cache=True
    )
    mol3 = stk.BuildingBlock(
        mol=join(bb_dir, 'neutral.mol'),
        functional_groups=['amine'],
        use_cache=True
    )
    assert mol0 is not mol2 and mol1 is not mol2
    assert mol2 is mol3

    mol4 = stk.BuildingBlock(
        mol=join(bb_dir, 'neutral.mol'),
        functional_groups=['aldehyde'],
        use_cache=True
    )
    assert mol3 is not mol4

    # Make sure that charged molecules are handled correctly.
    mol5 = stk.BuildingBlock(
        mol=join(bb_dir, 'negative_carbon.mol'),
        functional_groups=['amine'],
        use_cache=True
    )
    assert mol5 is not mol0
    # Test that all values are initialized correctly.
    assert len(mol5.func_groups) == 2
    fg_types = stk.dedupe(fg.info.name for fg in mol5.func_groups)
    assert sum(1 for _ in fg_types) == 1
    assert len(mol5.atoms) == 13
    assert len(mol5.bonds) == 12

    atom_count = {
        (stk.C, 0): 2,
        (stk.C, -1): 1,
        (stk.N, 0): 2,
        (stk.H, 0): 8,
    }
    assert atom_count == Counter(
        (a.__class__, a.charge) for a in mol5.atoms
    )

    expected_bonds = {
        frozenset({stk.N, stk.C}): 2,
        frozenset({stk.C}): 2,
        frozenset({stk.H, stk.N}): 4,
        frozenset({stk.H, stk.C}): 4
    }
    assert expected_bonds == Counter(
        frozenset({b.atom1.__class__, b.atom2.__class__})
        for b in mol5.bonds
    )

    mol6 = stk.BuildingBlock(
        mol=join(bb_dir, 'negative_nitrogen.mol'),
        functional_groups=['amine'],
        use_cache=True
    )
    assert mol6 is not mol5 and mol6 is not mol0
    # Test that all values are initialized correctly.
    assert len(mol6.func_groups) == 1
    fg_types = stk.dedupe(fg.info.name for fg in mol6.func_groups)
    assert sum(1 for _ in fg_types) == 1
    assert len(mol6.atoms) == 13
    assert len(mol6.bonds) == 12

    atom_count = {
        (stk.C, 0): 3,
        (stk.N, 0): 1,
        (stk.N, -1): 1,
        (stk.H, 0): 8,
    }
    assert atom_count == Counter(
        (a.__class__, a.charge) for a in mol6.atoms
    )

    expected_bonds = {
        frozenset({stk.N, stk.C}): 2,
        frozenset({stk.C}): 2,
        frozenset({stk.H, stk.N}): 2,
        frozenset({stk.H, stk.C}): 6
    }
    assert expected_bonds == Counter(
        frozenset({b.atom1.__class__, b.atom2.__class__})
        for b in mol6.bonds
    )


def test_init_pdb(bb_dir):
    mol0 = stk.BuildingBlock(join(bb_dir, 'neutral.pdb'), ['amine'])
    # Test that all values are initialized correctly.
    assert len(mol0.func_groups) == 2
    fg_types = stk.dedupe(fg.info.name for fg in mol0.func_groups)
    assert sum(1 for _ in fg_types) == 1
    assert len(mol0.atoms) == 15
    assert len(mol0.bonds) == 14

    atom_count = {
        (stk.H, 0): 10,
        (stk.N, 0): 2,
        (stk.C, 0): 3
    }
    assert atom_count == Counter(
        (a.__class__, a.charge) for a in mol0.atoms
    )

    expected_bonds = {
        frozenset({stk.N, stk.C}): 2,
        frozenset({stk.C}): 2,
        frozenset({stk.H, stk.N}): 4,
        frozenset({stk.H, stk.C}): 6
    }
    assert expected_bonds == Counter(
        frozenset({b.atom1.__class__, b.atom2.__class__})
        for b in mol0.bonds
    )

    # Test that caching is working properly.
    mol1 = stk.BuildingBlock(join(bb_dir, 'neutral.pdb'), ['amine'])
    assert mol0 is not mol1

    mol2 = stk.BuildingBlock(
        mol=join(bb_dir, 'neutral.pdb'),
        functional_groups=['amine'],
        use_cache=True
    )
    mol3 = stk.BuildingBlock(
        mol=join(bb_dir, 'neutral.pdb'),
        functional_groups=['amine'],
        use_cache=True
    )
    assert mol0 is not mol2 and mol1 is not mol2
    assert mol2 is mol3

    mol4 = stk.BuildingBlock(
        mol=join(bb_dir, 'neutral.pdb'),
        functional_groups=['aldehyde'],
        use_cache=True
    )
    assert mol3 is not mol4

    # Make sure that charged molecules are handled correctly.
    mol5 = stk.BuildingBlock(
        mol=join(bb_dir, 'negative_carbon.pdb'),
        functional_groups=['amine'],
        use_cache=True
    )
    assert mol5 is not mol0
    # Test that all values are initialized correctly.
    assert len(mol5.func_groups) == 2
    fg_types = stk.dedupe(fg.info.name for fg in mol5.func_groups)
    assert sum(1 for _ in fg_types) == 1
    assert len(mol5.atoms) == 13
    assert len(mol5.bonds) == 12

    atom_count = {
        (stk.C, 0): 2,
        (stk.C, -1): 1,
        (stk.N, 0): 2,
        (stk.H, 0): 8,
    }
    assert atom_count == Counter(
        (a.__class__, a.charge) for a in mol5.atoms
    )

    expected_bonds = {
        frozenset({stk.N, stk.C}): 2,
        frozenset({stk.C}): 2,
        frozenset({stk.H, stk.N}): 4,
        frozenset({stk.H, stk.C}): 4
    }
    assert expected_bonds == Counter(
        frozenset({b.atom1.__class__, b.atom2.__class__})
        for b in mol5.bonds
    )

    mol6 = stk.BuildingBlock(
        mol=join(bb_dir, 'negative_nitrogen.pdb'),
        functional_groups=['amine'],
        use_cache=True
    )
    assert mol6 is not mol5 and mol6 is not mol0
    # Test that all values are initialized correctly.
    assert len(mol6.func_groups) == 1
    fg_types = stk.dedupe(fg.info.name for fg in mol6.func_groups)
    assert sum(1 for _ in fg_types) == 1
    assert len(mol6.atoms) == 13
    assert len(mol6.bonds) == 12

    atom_count = {
        (stk.C, 0): 3,
        (stk.N, 0): 1,
        (stk.N, -1): 1,
        (stk.H, 0): 8,
    }
    assert atom_count == Counter(
        (a.__class__, a.charge) for a in mol6.atoms
    )

    expected_bonds = {
        frozenset({stk.N, stk.C}): 2,
        frozenset({stk.C}): 2,
        frozenset({stk.H, stk.N}): 2,
        frozenset({stk.H, stk.C}): 6
    }
    assert expected_bonds == Counter(
        frozenset({b.atom1.__class__, b.atom2.__class__})
        for b in mol6.bonds
    )


def test_init_from_random_file(bb_dir):
    mol0 = stk.BuildingBlock.init_from_random_file(
        file_glob=join(bb_dir, 'neutral.mol'),
        functional_groups=['amine']
    )
    # Test that all values are initialized correctly.
    assert len(mol0.func_groups) == 2
    fg_types = stk.dedupe(fg.info.name for fg in mol0.func_groups)
    assert sum(1 for _ in fg_types) == 1
    assert len(mol0.atoms) == 15
    assert len(mol0.bonds) == 14

    atom_count = {
        (stk.H, 0): 10,
        (stk.N, 0): 2,
        (stk.C, 0): 3
    }
    assert atom_count == Counter(
        (a.__class__, a.charge) for a in mol0.atoms
    )

    expected_bonds = {
        frozenset({stk.N, stk.C}): 2,
        frozenset({stk.C}): 2,
        frozenset({stk.H, stk.N}): 4,
        frozenset({stk.H, stk.C}): 6
    }
    assert expected_bonds == Counter(
        frozenset({b.atom1.__class__, b.atom2.__class__})
        for b in mol0.bonds
    )

    # Test that caching is working properly.
    mol1 = stk.BuildingBlock.init_from_random_file(
        file_glob=join(bb_dir, 'neutral.mol'),
        functional_groups=['amine']
    )
    assert mol0 is not mol1

    mol2 = stk.BuildingBlock.init_from_random_file(
        file_glob=join(bb_dir, 'neutral.mol'),
        functional_groups=['amine'],
        use_cache=True
    )
    mol3 = stk.BuildingBlock.init_from_random_file(
        file_glob=join(bb_dir, 'neutral.mol'),
        functional_groups=['amine'],
        use_cache=True
    )
    assert mol0 is not mol2 and mol1 is not mol2
    assert mol2 is mol3

    mol4 = stk.BuildingBlock.init_from_random_file(
        file_glob=join(bb_dir, 'neutral.mol'),
        functional_groups=['aldehyde'],
        use_cache=True
    )
    assert mol3 is not mol4

    # Make sure that charged molecules are handled correctly.
    mol5 = stk.BuildingBlock.init_from_random_file(
        file_glob=join(bb_dir, 'negative_carbon.mol'),
        functional_groups=['amine'],
        use_cache=True
    )
    assert mol5 is not mol0
    # Test that all values are initialized correctly.
    assert len(mol5.func_groups) == 2
    fg_types = stk.dedupe(fg.info.name for fg in mol5.func_groups)
    assert sum(1 for _ in fg_types) == 1
    assert len(mol5.atoms) == 13
    assert len(mol5.bonds) == 12

    atom_count = {
        (stk.C, 0): 2,
        (stk.C, -1): 1,
        (stk.N, 0): 2,
        (stk.H, 0): 8,
    }
    assert atom_count == Counter(
        (a.__class__, a.charge) for a in mol5.atoms
    )

    expected_bonds = {
        frozenset({stk.N, stk.C}): 2,
        frozenset({stk.C}): 2,
        frozenset({stk.H, stk.N}): 4,
        frozenset({stk.H, stk.C}): 4
    }
    assert expected_bonds == Counter(
        frozenset({b.atom1.__class__, b.atom2.__class__})
        for b in mol5.bonds
    )

    mol6 = stk.BuildingBlock.init_from_random_file(
        file_glob=join(bb_dir, 'negative_nitrogen.mol'),
        functional_groups=['amine'],
        use_cache=True
    )
    assert mol6 is not mol5 and mol6 is not mol0
    # Test that all values are initialized correctly.
    assert len(mol6.func_groups) == 1
    fg_types = stk.dedupe(fg.info.name for fg in mol6.func_groups)
    assert sum(1 for _ in fg_types) == 1
    assert len(mol6.atoms) == 13
    assert len(mol6.bonds) == 12

    atom_count = {
        (stk.C, 0): 3,
        (stk.N, 0): 1,
        (stk.N, -1): 1,
        (stk.H, 0): 8,
    }
    assert atom_count == Counter(
        (a.__class__, a.charge) for a in mol6.atoms
    )

    expected_bonds = {
        frozenset({stk.N, stk.C}): 2,
        frozenset({stk.C}): 2,
        frozenset({stk.H, stk.N}): 2,
        frozenset({stk.H, stk.C}): 6
    }
    assert expected_bonds == Counter(
        frozenset({b.atom1.__class__, b.atom2.__class__})
        for b in mol6.bonds
    )


def test_init_from_smiles():
    mol0 = stk.BuildingBlock.init_from_smiles('NCCCN', ['amine'])
    # Test that all values are initialized correctly.
    assert len(mol0.func_groups) == 2
    fg_types = stk.dedupe(fg.info.name for fg in mol0.func_groups)
    assert sum(1 for _ in fg_types) == 1
    assert len(mol0.atoms) == 15
    assert len(mol0.bonds) == 14

    atom_count = {
        (stk.H, 0): 10,
        (stk.N, 0): 2,
        (stk.C, 0): 3
    }
    assert atom_count == Counter(
        (a.__class__, a.charge) for a in mol0.atoms
    )

    expected_bonds = {
        frozenset({stk.N, stk.C}): 2,
        frozenset({stk.C}): 2,
        frozenset({stk.H, stk.N}): 4,
        frozenset({stk.H, stk.C}): 6
    }
    assert expected_bonds == Counter(
        frozenset({b.atom1.__class__, b.atom2.__class__})
        for b in mol0.bonds
    )

    # Test that caching is working properly.
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

    # Make sure that charged molecules are handled correctly.
    mol5 = stk.BuildingBlock.init_from_smiles(
        smiles='NC[C-]CN',
        functional_groups=['amine'],
        use_cache=True
    )
    assert mol5 is not mol0
    # Test that all values are initialized correctly.
    assert len(mol5.func_groups) == 2
    fg_types = stk.dedupe(fg.info.name for fg in mol5.func_groups)
    assert sum(1 for _ in fg_types) == 1
    assert len(mol5.atoms) == 13
    assert len(mol5.bonds) == 12

    atom_count = {
        (stk.C, 0): 2,
        (stk.C, -1): 1,
        (stk.N, 0): 2,
        (stk.H, 0): 8,
    }
    assert atom_count == Counter(
        (a.__class__, a.charge) for a in mol5.atoms
    )

    expected_bonds = {
        frozenset({stk.N, stk.C}): 2,
        frozenset({stk.C}): 2,
        frozenset({stk.H, stk.N}): 4,
        frozenset({stk.H, stk.C}): 4,
    }
    assert expected_bonds == Counter(
        frozenset({b.atom1.__class__, b.atom2.__class__})
        for b in mol5.bonds
    )

    mol6 = stk.BuildingBlock.init_from_smiles(
        smiles='[N-]CCCN',
        functional_groups=['amine'],
        use_cache=True
    )
    assert mol6 is not mol5 and mol6 is not mol0
    # Test that all values are initialized correctly.
    assert len(mol6.func_groups) == 1
    fg_types = stk.dedupe(fg.info.name for fg in mol6.func_groups)
    assert sum(1 for _ in fg_types) == 1
    assert len(mol6.atoms) == 13
    assert len(mol6.bonds) == 12

    atom_count = {
        (stk.C, 0): 3,
        (stk.N, 0): 1,
        (stk.N, -1): 1,
        (stk.H, 0): 8,
    }
    assert atom_count == Counter(
        (a.__class__, a.charge) for a in mol6.atoms
    )

    expected_bonds = {
        frozenset({stk.N, stk.C}): 2,
        frozenset({stk.C}): 2,
        frozenset({stk.H, stk.N}): 2,
        frozenset({stk.H, stk.C}): 6
    }
    assert expected_bonds == Counter(
        frozenset({b.atom1.__class__, b.atom2.__class__})
        for b in mol6.bonds
    )


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
    centroids = tmp_aldehyde3.get_bonder_centroids(fg_ids=[0])
    for i, centroid in enumerate(centroids):
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


def test_get_bonder_plane_normal(tmp_amine4):
    bonder_ids = list(tmp_amine4.get_bonder_ids())
    other_ids = [
        id_ for id_ in range(len(tmp_amine4.atoms))
        if id_ not in bonder_ids
    ]
    coords = tmp_amine4.get_position_matrix()
    coords[bonder_ids[0]] = [0, -1, 0]
    coords[bonder_ids[1]] = [-400, 0, 0]
    coords[bonder_ids[2]] = [0, 1, 1]
    coords[bonder_ids[3]] = [0, 1, -1]
    # Set the centroid of the molecule so that the plane normal
    # has a positive direction.
    coords[other_ids, 0] = 10
    coords[other_ids, 1] = 10
    tmp_amine4.set_position_matrix(coords)
    assert np.allclose(
        a=tmp_amine4.get_bonder_plane_normal(fg_ids=[0, 2, 3]),
        b=[1, 0, 0],
        atol=1e-6
    )

    # When using all fgs the plane should be at a 90 degree rotation.
    assert np.allcose(
        a=tmp_amine4.get_bonder_plane_normal(),
        b=[0, 1, 0],
        atol=1e-6
    )


def test_get_bonder_distances(tmp_amine4):
    # Place all bonders on a line.
    coords = tmp_amine4.get_position_matrix()
    for bonder_id in tmp_amine4.get_bonder_ids():
        coords[bonder_id] = [bonder_id, 0, 0]
    tmp_amine4.set_position_matrix(coords)

    # Test default behaviour.
    distances = tmp_amine4.get_bonder_distances()
    for i, (fg1, fg2, distance) in enumerate(distances):
        coord1 = tmp_amine4.func_groups[fg1].bonder_ids[0]
        coord2 = tmp_amine4.func_groups[fg2].bonder_ids[0]
        assert abs(distance - abs(coord1 - coord2)) < 1e-6
    assert i == 5

    # Test explicilty setting fg_ids.
    distances = tmp_amine4.get_bonder_distances(fg_ids=[0, 2, 3])
    for i, (fg1, fg2, distance) in enumerate(distances):
        coord1 = tmp_amine4.func_groups[fg1].bonder_ids[0]
        coord2 = tmp_amine4.func_groups[fg2].bonder_ids[0]
        assert abs(distance - abs(coord1 - coord2)) < 1e-6
    assert i == 2


def test_get_bonder_direction_vectors(tmp_amine4):
    pos_mat = tmp_amine4.get_position_matrix()
    # Set the coordinate of each bonder to the id of the fg.
    for fg in tmp_amine4.func_groups:
        for bonder in fg.bonder_ids:
            pos_mat[bonder] = [fg.id, fg.id, fg.id]
    tmp_amine4.set_position_matrix(pos_mat)

    dir_vectors = tmp_amine4.get_bonder_direction_vectors()
    for i, (id1, id2, v) in enumerate(dir_vectors):
        # Calculate the expected direction vector based on ids.
        d = stk.normalize_vector(np.array([id2]*3) - np.array([id1]*3))
        assert np.allclose(d, v, atol=1e-8)
    assert i == 5

    # Test explicitly setting fg_ids.
    dir_vectors = tmp_amine4.get_bonder_direction_vectors(
        fg_ids=[0, 3]
    )
    for i, (id1, id2, v) in enumerate(dir_vectors):
        # Calculate the expected direction vector based on ids.
        d = stk.normalize_vector(np.array([id2]*3) - np.array([id1]*3))
        assert np.allclose(d, v, atol=1e-8)
    assert i == 0


def test_get_centroid_centroid_direction_vector(amine4):
    c1 = amine4.get_centroid(atom_ids=amine4.get_bonder_ids())
    c2 = amine4.get_centroid()
    assert np.allclose(
        a=stk.normalize_vector(c2-c1),
        b=amine4.get_centroid_centroid_direction_vector(),
        atol=1e-8
    )

    # Test explicitly setting the fg_ids.
    fg_ids = [0, 2]
    c1 = amine4.get_centroid(atom_ids=amine4.get_bonder_ids(fg_ids))
    assert np.allclose(
        a=stk.normalize_vector(c2-c1),
        b=amine4.get_centroid_centroid_direction_vector(fg_ids=fg_ids),
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
    for a1, a2 in zip(tmp_amine2.atoms, mol2.atoms):
        assert vars(a1) == vars(a2)

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
