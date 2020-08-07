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
        np.allclose(i, j, atol=1e-4)
        for i, j in zip(
            periodic_case.cell,
            periodic_case.periodic_info.get_cell_matrix(),
        )
    ]))
