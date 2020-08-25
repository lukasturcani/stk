import numpy as np


def test_get_gamma(periodic_case):
    """
    Test :meth:`.PeriodicInfo.get_gamma`.

    Parameters
    ----------
    periodic_case : :class:`.CaseData`
        The periodic case to test.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert np.isclose(
        periodic_case.gamma,
        periodic_case.periodic_info.get_gamma(),
        atol=1e-6,
    )
