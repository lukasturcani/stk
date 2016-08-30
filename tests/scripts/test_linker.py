import os
import rdkit
import numpy as np

from ...classes import Linker, FGInfo
from .test_struct_unit import get_mol_file

lk_name = next(x for x in get_mol_file() 
                                if 'test_rot_amine.mol' in x)  
lk = Linker(lk_name)

def test_get_heavy_direction_vector():

    assert np.array_equal(lk.get_heavy_direction_vector(), 
                          np.array([1, 0, 0]))
                          
def test_get_heavy_theta():
    
    # Theta between the vectors defining the axes should be 90 degrees,
    # except for the x-axis which should be 0 degrees.
    assert lk.get_heavy_theta(np.array([1, 0, 0])) == 0
    assert lk.get_heavy_theta(np.array([0, 1, 0])) == np.pi/2
    assert lk.get_heavy_theta(np.array([0, 0, 1])) == np.pi/2
    


    
    
    
    
    
    
    