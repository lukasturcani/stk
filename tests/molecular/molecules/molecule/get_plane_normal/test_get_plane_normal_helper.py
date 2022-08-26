import numpy as np
import pytest

from .utilities import get_plane_normal


@pytest.mark.parametrize(
    argnames=("position_matrix", "atom_ids", "normal"),
    argvalues=(
        (
            np.array(
                [
                    [1.0, 2.0, 0.0],
                    [3.0, 4.0, 0.0],
                    [123.0, 423, 0.0],
                ]
            ),
            (0, 1, 2),
            np.array([0.0, 0.0, 1.0]),
        ),
        (
            np.array(
                [
                    [1.0, 2.0, 0.0],
                    [3.0, 4.0, 0.0],
                    [123.0, 423, 105.0],
                ]
            ),
            (0, 1),
            np.array([0.0, 0.0, 1.0]),
        ),
        (
            np.array(
                [
                    [0.0, 2.0, 3.0],
                    [0.0, 4.0, 1.0],
                    [0.0, 423, 12.0],
                ]
            ),
            (0, 1, 2),
            np.array([1.0, 0.0, 0.0]),
        ),
    ),
)
def test_get_plane_normal_helper(position_matrix, atom_ids, normal):
    """
    Test :func:`.get_plane_normal`.

    Parameters
    ----------
    position_matrix : :class:`numpy.ndarray`
        The position matrix to test.

    atom_ids : :class:`tuple` of :class:`int`
        The ids of atoms to include in the test.

    normal : :class:`numpy.ndarray`
        The correct plane normal.

    Returns
    -------
    None : :class:`NoneType`

    """

    result = get_plane_normal(position_matrix, atom_ids)
    # The normal may be parallel or anti-parallel.
    assert abs(abs(result @ normal) - 1) < 1e-32
