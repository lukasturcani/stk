import numpy as np
import os
import pickle

from ...convenience_functions import normalize_vector, vector_theta

mol_file = os.path.join('data', 'linker', 'linker_test_obj')
with open(mol_file, 'rb') as dump_file:
    lk = pickle.load(dump_file)

def test_set_heavy_mol_orientation():
    """
    Tests `set_heavy_mol_orientation`.
    
    """
    
    # Set it and check it.       
    lk.set_heavy_mol_orientation([1,2,3])
    assert np.allclose(next(lk.heavy_direction_vectors()),
                       normalize_vector([1,2,3]),
                       atol=1e-8)
                       
def test_minimize_theta():
    """
    Tests `minimize_theta`.
        
    """

    with open(mol_file, 'rb') as dump_file:
        lk = pickle.load(dump_file)      
    
    vector = [1,1,1]
    
    dir_vector = next(lk.heavy_direction_vectors())
    og_theta = vector_theta(vector, dir_vector)  
    lk.minimize_theta(vector, [0,0,1])
    dir_vector = next(lk.heavy_direction_vectors())
    new_theta = vector_theta(vector, dir_vector)
    assert new_theta < og_theta