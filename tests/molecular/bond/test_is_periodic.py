def test_is_periodic(case_data):
    """
    Test :meth:`.Bond.is_periodic`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the bond to test and its true periodicity.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_is_periodic(case_data.bond, case_data.periodicity)


def _test_is_periodic(bond, periodicity):
    """
    Test :meth:`.Bond.is_periodic`.

    Parameters
    ----------
    bond : :class:`.Bond`
        The bond to test.

    periodicity : :class:`tuple` of :class:`int`
        The true periodicity of `bond`.

    Returns
    -------
    None : :class:`NoneType`

    """

    return bond.is_periodic() == any(p != 0 for p in periodicity)
