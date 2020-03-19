def test_get_atom(case_data):
    """
    Test :meth:`.AtomInfo.get_atom`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the atom info to test and the correct atom.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert case_data.atom_info.get_atom() is case_data.atom
