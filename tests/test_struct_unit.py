import os
import numpy as np
import rdkit.Chem.AllChem as rdkit
import stk

if not os.path.exists('struct_unit_tests_output'):
    os.mkdir('struct_unit_tests_output')


def test_init(mol2):
    assert len(mol2.functional_group_atoms()) == 3
    assert mol2.func_grp.name == 'amine'


def test_all_bonder_distances(mol2):
    for i, (id1, id2, d) in enumerate(mol2.all_bonder_distances()):
        assert type(id1) == int
        assert type(id2) == int
        assert type(d) == float
    assert i == 2


def test_bonder_centroids(mol2):
    for i, centroid in enumerate(mol2.bonder_centroids()):
        assert len(centroid) == 3
        assert type(centroid) == np.ndarray
    assert i == 2


def test_bonder_centroid(mol2):
    centroid = mol2.bonder_centroid()
    assert len(centroid) == 3
    assert type(centroid) == np.ndarray


def test_bonder_direction_vectors(mol2):
    for i, (id1, id2, v) in enumerate(mol2.bonder_direction_vectors()):
        assert type(id1) == int
        assert type(id2) == int
        assert type(v) == np.ndarray
    assert i == 2


def test_bonder_position_matrix(mol2):
    m = mol2.bonder_position_matrix()
    assert len(mol2.bonder_ids) == m.shape[1]


def test_centroid_centroid_dir_vector(mol2):
    c1 = mol2.bonder_centroid()
    c2 = mol2.centroid()
    assert np.allclose(stk.normalize_vector(c2-c1),
                       mol2.centroid_centroid_dir_vector(),
                       atol=1e-8)


def test_core(mol2):
    for atom in mol2.core().GetAtoms():
        assert atom.GetAtomicNum() != 1
        assert not atom.HasProp('fg')


def test_functional_group_atoms(mol2):
        func_grp_mol = rdkit.MolFromSmarts(mol2.func_grp.fg_smarts)
        assert (mol2.mol.GetSubstructMatches(func_grp_mol) ==
                mol2.functional_group_atoms())


def test_is_core_atom(mol2):
    for atom in mol2.mol.GetAtoms():
        core = (False if atom.HasProp('fg') or atom.GetAtomicNum() == 1
                else True)
        assert core is mol2.is_core_atom(atom.GetIdx())


def test_json_init(mol):
    path = os.path.join('struct_unit_tests_output', 'mol.json')
    mol.dump(path)
    mol2 = stk.Molecule.load(path, stk.Molecule.from_dict)

    assert isinstance(mol.file, str)
    assert mol2.optimized
    assert mol2.bonder_ids == mol.bonder_ids
    assert mol2.energy.__class__.__name__ == 'Energy'
    assert mol2.func_grp.name == 'amine'
    assert mol is not mol2
    assert mol.atom_props == mol.atom_props


def test_caching():
    # This test is in a try block because pytest runs with
    # CACHE_SETTINGS turned off. If this test fails, the finally
    # clause ensures that the cache remains off so other tests are
    # not affected by unexpected cache problems.
    try:
        stk.CACHE_SETTINGS['ON'] = True
        mol = stk.StructUnit.smiles_init('NC1CC(N)CC(N)C1', 'amine')
        mol2 = stk.StructUnit.smiles_init('NC1CC(N)CC(N)C1', 'amine')
        assert mol is mol2

        mol3 = stk.StructUnit.smiles_init('NC1CC(N)CC(N)C1', 'aldehyde')
        assert mol3 is not mol

        stk.CACHE_SETTINGS['ON'] = False
        mol4 = stk.StructUnit.smiles_init('NC1CC(N)CC(N)C1', 'amine')
        stk.CACHE_SETTINGS['ON'] = True

        assert mol is not mol4

    except Exception:
        raise

    finally:
        stk.CACHE_SETTINGS['ON'] = False


def test_set_bonder_centroid(tmp_mol):
    tmp_mol.set_bonder_centroid([1, 2, 3])
    assert np.allclose(tmp_mol.bonder_centroid(), [1, 2, 3], atol=1e-8)


def test_untag_atoms(tmp_mol):
    assert any(a.HasProp('fg') for a in tmp_mol.mol.GetAtoms())
    tmp_mol.untag_atoms()
    assert all(not a.HasProp('fg') for a in tmp_mol.mol.GetAtoms())
