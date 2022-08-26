import numpy as np
import pytest

from .utilities import get_maximum_diameter


@pytest.mark.parametrize(
    argnames=("position_matrix", "atom_ids", "maximum_diameter"),
    argvalues=(
        (
            np.array(
                [
                    [1.0, 2.0, 3.0],
                    [1.0, 2.0, 4.0],
                ]
            ),
            (0,),
            0,
        ),
        (
            np.array(
                [
                    [1.0, 0.0, 0.0],
                    [-1.0, 0.0, 0.0],
                    [9.0, 0.0, 0.0],
                ]
            ),
            (0, 1, 2),
            10,
        ),
        (
            np.array(
                [
                    [1.0, 0.0, 0.0],
                    [-1.0, 0.0, 0.0],
                    [9.0, 0.0, 0.0],
                ]
            ),
            (0, 2),
            8,
        ),
        (
            np.array(
                [
                    [0.0, 0.0, 0.0],
                    [3.0, 4.0, 0.0],
                ]
            ),
            (0, 1),
            5,
        ),
    ),
)
def test_get_direction_helper(
    position_matrix,
    atom_ids,
    maximum_diameter,
):
    """
    Test :func:`.get_maximum_diameter`.

    Parameters
    ----------
    position_matrix : :class:`numpy.ndarray`
        The position matrix to test.

    atom_ids : :class:`tuple` of :class:`int`
        The ids of atoms included in the test.

    maximum_diameter : :class:`float`
        The correct maximum diameter.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert np.allclose(
        a=get_maximum_diameter(position_matrix, atom_ids),
        b=maximum_diameter,
        atol=1e-32,
    )
