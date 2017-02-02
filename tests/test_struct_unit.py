from os.path import join
import numpy as np
import itertools as it
from scipy.spatial.distance import euclidean
import rdkit.Chem as chem

from ..molecular import StructUnit
from ..convenience_tools import normalize_vector

data_dir = join('data', 'struct_unit', 'amine.mol')
mol = StructUnit(data_dir)
conf = mol.mol.GetConformer()

def test_init():
    f = join('data', 'struct_unit', 'amine2.mol2')
    mol2 = StructUnit(f)
    mol3 = StructUnit(f, 'aldehyde')
    print(mol2.func_grp.name)
    assert len(mol2.functional_group_atoms()) == 3
    assert len(mol3.functional_group_atoms()) == 2

def test_all_bonder_distances():
    conf_ds = []
    for a1, a2 in it.combinations(mol.bonder_ids, 2):
        a1_coord = conf.GetAtomPosition(a1)
        a2_coord = conf.GetAtomPosition(a2)
        conf_ds.append(euclidean(a1_coord, a2_coord))

    assert len(conf_ds) == sum(1 for _ in mol.all_bonder_distances())
    assert np.allclose(list(x for *_, x in mol.all_bonder_distances()), 
                       conf_ds, atol=1e-8)

def test_bonder_centroid():
    position = np.array([0.,0.,0.])
    for id_ in mol.bonder_ids:
        position += np.array([*conf.GetAtomPosition(id_)])
        
    assert np.allclose((position/len(mol.bonder_ids)), 
                       mol.bonder_centroid(), atol=1e-8)
                       
def test_bonder_direction_vectors():
    vs = []
    for atom1_id, atom2_id in it.combinations(mol.bonder_ids, 2):
        p1 = np.array([*conf.GetAtomPosition(atom1_id)])
        p2 = np.array([*conf.GetAtomPosition(atom2_id)])
        
        vs.append(normalize_vector(p1-p2))
        
    assert len(vs) == sum(1 for _ in mol.bonder_direction_vectors())
    assert np.allclose(
        list(x for *_, x in mol.bonder_direction_vectors()), vs,
            atol=1e-8)
            
def test_bonder_position_matrix():
        pos_array = np.array([])

        for atom_id in mol.bonder_ids:
            pos_vect = np.array([*conf.GetAtomPosition(atom_id)])
            pos_array = np.append(pos_array, pos_vect)

        m = np.matrix(pos_array.reshape(-1,3).T)
    
        assert np.allclose(m, mol.bonder_position_matrix(), atol=1e-8)
        
def test_centroid_centroid_dir_vector():
    c1 = mol.bonder_centroid()
    c2 = mol.centroid()
    assert np.allclose(normalize_vector(c2-c1), 
                       mol.centroid_centroid_dir_vector(),
                       atol=1e-8)
                       
def test_functional_group_atoms():  
        func_grp_mol = chem.MolFromSmarts(mol.func_grp.fg_smarts)
        assert (mol.mol.GetSubstructMatches(func_grp_mol)  == 
                mol.functional_group_atoms())

def test_set_bonder_centroid():
    og = mol.bonder_centroid()
    mol.set_bonder_centroid([1,2,3])
    assert np.allclose(mol.bonder_centroid(), [1,2,3], atol=1e-8)
    mol.set_bonder_centroid(og)

def test_untag_atoms():
    f = join('data', 'struct_unit', 'amine2.mol2')
    mol2 = StructUnit(f)
    mol2.untag_atoms()
    for atom in mol2.mol.GetAtoms():
        assert not atom.HasProp('fg')