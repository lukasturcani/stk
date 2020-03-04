import numpy as np


def test_get_position(test_case):
    _test_get_position(test_case.vertex, test_case.position)


def _test_get_position(vertex, position):
    assert np.all(np.equal(vertex.get_position(), position))
