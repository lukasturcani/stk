import numpy as np
import stk


def test_set_orientation2(tmp_aldehyde3):
    tmp_aldehyde3.set_orientation2([1, 2, 3], 0)
    assert np.allclose(tmp_aldehyde3.bonder_plane_normal(0),
                       stk.normalize_vector([1, 2, 3]), atol=1e-4)
