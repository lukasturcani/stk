import numpy as np


def test_get_alpha(periodic_case):
    """
    Test :meth:`.PeriodicInfo.get_alpha`.

    Parameters
    ----------
    periodic_case : :class:`.CaseData`
        The periodic case to test.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert np.isclose(
        periodic_case.alpha,
        periodic_case.periodic_info.get_alpha(),
        atol=1e-6,
    )
