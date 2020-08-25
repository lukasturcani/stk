import numpy as np


def test_get_a(periodic_case):
    """
    Test :meth:`.PeriodicInfo.get_a`.

    Parameters
    ----------
    periodic_case : :class:`.CaseData`
        The periodic case to test.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert np.isclose(
        periodic_case.a,
        periodic_case.periodic_info.get_a(),
        atol=1e-6,
    )
