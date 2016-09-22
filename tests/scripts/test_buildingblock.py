import numpy as np
import pickle
import os
from rdkit.Geometry.rdGeometry import Point3D

from MMEA.convenience_functions import normalize_vector

# These tests use a flat building block so that verification of the 
# plane calculations are easy.

bb_file = os.path.join('data', 'building_block', 
                       'building_block_test_obj')
with open(bb_file, 'rb') as dump_file:
    bb = pickle.load(dump_file)

def test_heavy_orientation_functions():
    """
    Tests functions used for setting the orientation.
    
    Tested functions:
        > heavy_plane_normal
        > set_heavy_mol_orientation
    
    """
    # This function does modifications so load a new copy of bb.
    with open(bb_file, 'rb') as dump_file:
        bb = pickle.load(dump_file)


    # Test object has the 3 heavy atoms on a flat plane along the xy
    # axis.
    assert np.allclose(bb.heavy_plane_normal(), [0,0,-1])
    assert not np.allclose(bb.heavy_plane_normal(), [0,0,1])
    # The normal should be pointing toward the general direction of the
    # centroid. In the test it is expected to be negative which means
    # the z position of the molecular centroid is lower than the z 
    # position of the heavy atom centroid.
    assert bb.centroid('heavy')[-1] < bb.heavy_atom_centroid()[-1]
    
    # Changing the position of the Carbon atom to be large and positive
    # should flip the direction of the heavy_plane_normal.
    conf = bb.heavy_mol.GetConformer()
    new_pos = Point3D(0,0,100)
    for atom in bb.heavy_mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            conf.SetAtomPosition(atom.GetIdx(), new_pos)
            break
    
    assert not np.allclose(bb.heavy_plane_normal(), [0,0,-1])
    assert np.allclose(bb.heavy_plane_normal(), [0,0,1])
    
    # Set it to something else and check that it worked.
    bb.set_heavy_mol_orientation([1,0,0])
    assert np.allclose(bb.heavy_plane_normal(), [1,0,0])
  
    
def test_heavy_plane():
    """
    Tests `heavy_plane`.
    
    """
    
    # The first 3 members of the plane array should be equal to the 
    # heavy plane normal.
    assert np.allclose(bb.heavy_plane()[:3], bb.heavy_plane_normal(), 
                       atol=1e-8)     

def test_centroid_centroid_dir_vector():
    """
    Tests `centroid_centroid_dir_vector`.
    
    """
    
    mol_centroid = bb.centroid('heavy')
    h_atom_centroid = bb.heavy_atom_centroid()
    
    assert np.array_equal(
        normalize_vector(mol_centroid - h_atom_centroid),
                          bb.centroid_centroid_dir_vector())