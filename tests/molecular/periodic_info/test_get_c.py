import numpy as np


def test_get_c(periodic_case):
    """
    Test :meth:`.PeriodicInfo.get_c`.

    Parameters
    ----------
    periodic_case : :class:`.CaseData`
        The periodic case to test.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert np.isclose(
        periodic_case.c,
        periodic_case.periodic_info.get_c(),
        atol=1e-6,
    )
