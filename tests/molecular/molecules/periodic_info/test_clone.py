import numpy as np


def test_clone(periodic_case):
    """
    Test :meth:`.PeriodicInfo.clone`.

    Parameters
    ----------
    periodic_case : :class:`.CaseData`
        The periodic case to test.

    Returns
    -------
    None : :class:`NoneType`

    """

    clone = periodic_case.periodic_info.clone()

    assert np.all(np.array([
        np.allclose(i, j, atol=1e-4)
        for i, j in zip(
            clone.get_cell_matrix(),
            periodic_case.periodic_info.get_cell_matrix()
        )
    ]))
