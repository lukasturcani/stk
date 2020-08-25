import numpy as np


def test_get_b(periodic_case):
    """
    Test :meth:`.PeriodicInfo.get_b`.

    Parameters
    ----------
    periodic_case : :class:`.CaseData`
        The periodic case to test.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert np.isclose(
        periodic_case.b,
        periodic_case.periodic_info.get_b(),
        atol=1e-6,
    )
