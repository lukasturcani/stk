import numpy as np


def test_get_cell_matrix(periodic_case):
    """
    Test :meth:`.PeriodicInfo.get_cell_matrix`.

    Parameters
    ----------
    periodic_case : :class:`.CaseData`
        The periodic case to test.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert np.all(np.array([
        np.allclose(i, j, atol=1e-6)
        for i, j in zip(
            (
                periodic_case.x_vector,
                periodic_case.y_vector,
                periodic_case.z_vector
            ),
            periodic_case.periodic_info.get_cell_matrix(),
        )
    ]))
