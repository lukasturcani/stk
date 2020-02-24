import numpy as np
import stk
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
def test_with_scale(start, scale, target):
    edge = stk.Edge(
        id=0,
        vertex1=stk.Vertex(0, [0, 0, 0]),
        vertex2=stk.Vertex(1, [10, 0, 0]),
        position=start,
    )
    # Test immutability.
    clone = edge.clone()
    _test_with_scale(edge, scale, target)
    is_clone(edge, clone)


def _test_with_scale(edge, scale, target):
    clone = edge.with_scale(scale)
    assert np.allclose(
        a=clone.get_position(),
        b=target,
        atol=1e-16,
    )
