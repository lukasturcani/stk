import numpy as np
from os.path import join

from ..convenience_tools import normalize_vector
from ..molecular import StructUnit3

mol_file = join('data', 'struct_unit3', 'amine.mol2')

def test_set_orientation2():
    mol = StructUnit3(mol_file)
    mol.set_orientation2([1,2,3])
    assert np.allclose(mol.bonder_plane_normal(), 
                           normalize_vector([1,2,3]), atol=1e-8)