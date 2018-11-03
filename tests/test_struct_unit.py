import os
import numpy as np
import rdkit.Chem.AllChem as rdkit
import stk

if not os.path.exists('struct_unit_tests_output'):
    os.mkdir('struct_unit_tests_output')

mol = stk.StructUnit.smiles_init('NC1CC(N)CC(N)C1', 'amine')
conf = mol.mol.GetConformer()


def test_init():
    assert len(mol.functional_group_atoms()) == 3
    assert mol.func_grp.name == 'amine'


def test_all_bonder_distances():
    for i, (id1, id2, d) in enumerate(mol.all_bonder_distances()):
        assert type(id1) == int
        assert type(id2) == int
        assert type(d) == float
    assert i == 2


def test_bonder_centroids():
    for i, centroid in enumerate(mol.bonder_centroids()):
        assert len(centroid) == 3
        assert type(centroid) == np.ndarray
    assert i == 2


def test_bonder_centroid():
    centroid = mol.bonder_centroid()
    assert len(centroid) == 3
    assert type(centroid) == np.ndarray


def test_bonder_direction_vectors():
    for i, (id1, id2, v) in enumerate(mol.bonder_direction_vectors()):
        assert type(id1) == int
        assert type(id2) == int
        assert type(v) == np.ndarray
    assert i == 2


def test_bonder_position_matrix():
    m = mol.bonder_position_matrix()
    assert len(mol.bonder_ids) == m.shape[1]


def test_centroid_centroid_dir_vector():
    c1 = mol.bonder_centroid()
    c2 = mol.centroid()
    assert np.allclose(stk.normalize_vector(c2-c1),
                       mol.centroid_centroid_dir_vector(),
                       atol=1e-8)


def test_core():
    for atom in mol.core().GetAtoms():
        assert atom.GetAtomicNum() != 1
        assert not atom.HasProp('fg')


def test_functional_group_atoms():
        func_grp_mol = rdkit.MolFromSmarts(mol.func_grp.fg_smarts)
        assert (mol.mol.GetSubstructMatches(func_grp_mol) ==
                mol.functional_group_atoms())


def test_is_core_atom():
    for atom in mol.mol.GetAtoms():
        core = (False if atom.HasProp('fg') or atom.GetAtomicNum() == 1
                else True)
        assert core is mol.is_core_atom(atom.GetIdx())


def test_json_init():
    try:
        path = os.path.join('struct_unit_tests_output', 'mol.json')
        mol.dump(path)
        stk.CACHE_SETTINGS['ON'] = False
        mol2 = stk.Molecule.load(path, stk.Molecule.from_dict)
        stk.CACHE_SETTINGS['ON'] = True

        assert isinstance(mol.file, str)
        assert mol2.optimized
        assert mol2.bonder_ids == mol.bonder_ids
        assert mol2.energy.__class__.__name__ == 'Energy'
        assert mol2.func_grp.name == 'amine'
        assert mol is not mol2
        assert mol.atom_props == mol.atom_props

    except Exception:
        raise
    finally:
        stk.CACHE_SETTINGS['ON'] = True


def test_caching():
    try:
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
        stk.CACHE_SETTINGS['ON'] = True


def test_set_bonder_centroid():
    mol.set_bonder_centroid([1, 2, 3])
    assert np.allclose(mol.bonder_centroid(), [1, 2, 3], atol=1e-8)


def test_untag_atoms():
    try:
        stk.CACHE_SETTINGS['ON'] = False
        mol = stk.StructUnit.smiles_init('NC1CC(N)CC(N)C1', 'amine')
        stk.CACHE_SETTINGS['ON'] = True
        assert any(a.HasProp('fg') for a in mol.mol.GetAtoms())
        mol.untag_atoms()
        assert all(not a.HasProp('fg') for a in mol.mol.GetAtoms())
    except Exception:
        raise
    finally:
        stk.CACHE_SETTINGS['ON'] = True
