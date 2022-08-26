import numpy as np
import pytest

from ..utilities import get_centroid


@pytest.mark.parametrize(
    argnames=("position_matrix", "atom_ids", "centroid"),
    argvalues=(
        (
            np.array(
                [
                    [1.0, 2.0, 3.0],
                    [0.0, 0.0, 0.0],
                    [-1.0, -2.0, -3.0],
                ]
            ),
            (0, 1, 2),
            np.array([0.0, 0.0, 0.0]),
        ),
        (
            np.array(
                [
                    [1.0, 2.0, 3.0],
                    [0.0, 0.0, 0.0],
                    [-1.0, -2.0, -3.0],
                ]
            ),
            (0, 1),
            np.array([0.5, 1.0, 1.5]),
        ),
    ),
)
def test_get_centroid_helper(position_matrix, atom_ids, centroid):
    """
    Test :func:`.get_centroid`.

    Parameters
    ----------
    position_matrix : :class:`numpy.ndarray`
        A position matrix.

    atom_ids : :class:`tuple` of :class:`int`
        The ids of selected atoms.

    centroid : :class:`numpy.ndarray`
        The correct centroid.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert np.allclose(
        a=get_centroid(position_matrix, atom_ids),
        b=centroid,
        atol=1e-32,
    )
