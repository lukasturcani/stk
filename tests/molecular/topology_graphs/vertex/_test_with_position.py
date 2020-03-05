import numpy as np

from .utilities import is_clone


def test_with_position(case_data, position):
    _test_with_position(case_data.vertex, position)


def _test_with_position(vertex, position):
    # Test immutability.
    clone = vertex.clone()
    _test_with_position_0(vertex, position)
    is_clone(vertex, clone)


def _test_with_position_0(vertex, position):
    clone = vertex.with_position(position)
    assert np.all(np.equal(clone.get_position(), position))
