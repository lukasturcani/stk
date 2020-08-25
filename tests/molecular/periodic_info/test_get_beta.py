import numpy as np


def test_get_beta(periodic_case):
    """
    Test :meth:`.PeriodicInfo.get_beta`.

    Parameters
    ----------
    periodic_case : :class:`.CaseData`
        The periodic case to test.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert np.isclose(
        periodic_case.beta,
        periodic_case.periodic_info.get_beta(),
        atol=1e-6,
    )
