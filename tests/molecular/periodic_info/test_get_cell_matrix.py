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

    test = np.array(
        (
            periodic_case.vector_1,
            periodic_case.vector_2,
            periodic_case.vector_3,
        )
    )
    original = np.array(periodic_case.periodic_info.get_cell_matrix())

    assert np.all(np.allclose(test, original, atol=1e-6))
