def test_get_charge(case_data):
    """
    Test :meth:`.Atom.get_charge`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        The test case. Holds the atom to test and the expected result.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert case_data.atom.get_charge() == case_data.charge
