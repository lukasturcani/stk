import numpy as np
import stk


def test_set_orientation2(mol4):
    mol4.set_orientation2([1, 2, 3])
    assert np.allclose(mol4.bonder_plane_normal(),
                       stk.normalize_vector([1, 2, 3]), atol=1e-4)
