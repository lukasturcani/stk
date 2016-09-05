import numpy as np

from ...classes import Linker
from .test_struct_unit import get_mol_file
from ...convenience_functions import normalize_vector

lk_name = next(x for x in get_mol_file() 
                                if 'test_rot_amine.mol' in x)  
lk = Linker(lk_name)

def test_set_heavy_mol_orientation():
    lk.set_heavy_mol_orientation([1,2,3])
    assert np.allclose(next(lk.heavy_direction_vectors()),
                       normalize_vector([1,2,3]),
                       atol=1e-8)
                       
    