import numpy as np
from os.path import join

from ..convenience_tools import normalize_vector
from ..molecular import StructUnit2

mol_file = join('data', 'struct_unit2', 'amine.mol2')

def test_set_orientation2():
    mol = StructUnit2(mol_file)
    mol.set_orientation2([1,2,3])
    assert np.allclose(next(mol.bonder_direction_vectors())[-1], 
                           normalize_vector([1,2,3]), atol=1e-8)