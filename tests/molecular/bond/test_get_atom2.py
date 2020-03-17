def test_get_atom2(case_data):
    """
    Test :meth:`.Bond.get_atom2`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the bond to test and the correct *atom2*.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert case_data.bond.get_atom2() is case_data.atom2
