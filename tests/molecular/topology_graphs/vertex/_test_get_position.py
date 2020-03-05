import numpy as np


def test_get_position(case_data):
    _test_get_position(case_data.vertex, case_data.position)


def _test_get_position(vertex, position):
    assert np.allclose(vertex.get_position(), position, atol=1e-15)
