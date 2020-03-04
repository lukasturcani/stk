import numpy as np
import pytest

from .utilities import is_clone


@pytest.mark.parametrize(
    argnames=('start', 'scale', 'target'),
    argvalues=(
        (np.array([0, 0, 0]), 2, np.array([0, 0, 0])),
        (np.array([0, 0, 0]), 0, np.array([0, 0, 0])),
        (np.array([1, -1, 0]), 0, np.array([0, 0, 0])),
        (np.array([1, -1, 2]), 2, np.array([2, -2, 4])),
        (np.array([2, 4, -5]), -3, np.array([-6, -12, 15])),
        (
            np.array([1, -1, 2]),
            np.array([10, 2, -2]),
            np.array([10, -2, -4]),
        ),
    ),
)
def test_with_scale(test_case, start, scale, target):
    _test_with_scale(
        vertex=test_case.vertex,
        start=start,
        scale=scale,
        target=target,
    )


def _test_with_scale(vertex, start, scale, target):
    # Test immutability.
    clone = vertex.clone()
    _test_with_scale_0(vertex, start, scale, target)
    is_clone(vertex, clone)


def _test_with_scale_0(vertex, start, scale, target):
    clone = vertex.with_position(start).with_scale(scale)
    assert np.allclose(
        a=clone.get_position(),
        b=target,
        atol=1e-16,
    )
