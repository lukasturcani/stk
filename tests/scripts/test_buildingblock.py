import numpy as np

from ...classes import BuildingBlock
from .test_struct_unit import get_mol_file
from ...convenience_functions import normalize_vector

bb_name = next(x for x in get_mol_file() 
         if 'buildingblock_flat_amine_test.mol' in x)  
bb = BuildingBlock(bb_name)

def test_heavy_plane_normal():
    
    assert np.allclose(bb.heavy_plane_normal(), [0,0,1])
    
    bb.set_heavy_mol_orientation([1,0,0])
    
    assert np.allclose(bb.heavy_plane_normal(), [1,0,0])
    
    bb.set_heavy_mol_orientation([0,0,1])    
    
def test_heavy_plane():
        

    # The first 3 members of the plane array should be equal to the 
    # heavy plane normal.
    assert np.allclose(bb.heavy_plane()[:3], bb.heavy_plane_normal(), 
                       atol=1e-8)     

def test_set_heavy_mol_orientation():
    
    assert np.allclose(bb.heavy_plane_normal(), [0,0,1])
    
    bb.set_heavy_mol_orientation([1,1,1])
    
    assert np.allclose(bb.heavy_plane_normal(), 
                       normalize_vector([1,1,1]), atol=1e-3)
    
    bb.set_heavy_mol_orientation([2,3,10])
    
    assert np.allclose(bb.heavy_plane_normal(), 
                       normalize_vector([2,3,10]), atol=1e-3)
                       
    bb.set_heavy_mol_orientation([0,0,1])