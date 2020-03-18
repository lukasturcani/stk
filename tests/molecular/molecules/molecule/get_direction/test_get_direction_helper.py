import pytest
import numpy as np

from ..utilities import get_direction


@pytest.mark.parametrize(
    argnames=('position_matrix', 'atom_ids', 'direction'),
    argvalues=(
        (
            np.array([
                [1., 0., 0.],
                [2., 0., 0.],
                [3., 0., 0.],
            ]),
            (0, 1, 2),
            np.array([1., 0., 0.]),
        ),
        (
            np.array([
                [0., 1., 0.],
                [0., 2., 0.],
                [0., 3., 0.],
            ]),
            (0, 1, 2),
            np.array([0., 1., 0.]),
        ),
        (
            np.array([
                [1., 0., 0.],
                [2., 0., 0.],
                [0., 5., 0.],
            ]),
            (0, 1),
            np.array([1., 0., 0.]),
        ),
    ),
)
def test_get_direction_helper(position_matrix, atom_ids, direction):
    """
    Test the :func:`.get_direction`.

    """

    assert np.allclose(
        a=get_direction(position_matrix, atom_ids),
        b=direction,
        atol=1e-32,
    )
