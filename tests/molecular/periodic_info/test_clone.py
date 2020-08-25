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

    test = np.array(clone.get_cell_matrix())
    original = np.array(
        periodic_case.periodic_info.get_cell_matrix()
    )

    assert np.all(np.equal(test, original))
