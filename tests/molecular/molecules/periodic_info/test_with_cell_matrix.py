import numpy as np


def test_with_cell_matrix(periodic_case):
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

    new_info = periodic_case.periodic_info.with_cell_matrix(
        cell_matrix=periodic_case.cell
    )

    assert np.all(np.array([
        np.allclose(i, j, atol=1e-4)
        for i, j in zip(
            periodic_case.cell,
            new_info.get_cell_matrix(),
        )
    ]))
