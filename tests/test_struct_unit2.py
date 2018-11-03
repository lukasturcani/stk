import numpy as np
import stk

mol = stk.StructUnit2.smiles_init('NCCCN', 'amine')


def test_set_orientation2():
    mol.set_orientation2([1, 2, 3])
    assert np.allclose(next(mol.bonder_direction_vectors())[-1],
                       stk.normalize_vector([1, 2, 3]), atol=1e-8)
