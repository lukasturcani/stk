import numpy as np
import stk


def test_set_orientation2(tmp_amine2):
    tmp_amine2.set_orientation2([1, 2, 3], 0)
    vector = next(tmp_amine2.bonder_direction_vectors(0))[-1]
    assert np.allclose(vector,
                       stk.normalize_vector([1, 2, 3]),
                       atol=1e-8)
