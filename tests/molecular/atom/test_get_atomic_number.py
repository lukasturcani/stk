def test_get_atomic_number(case_data):
    """
    Test :meth:`.Atom.get_atomic_number`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case, containing the atom to test and its correct atomic
        number.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_get_atomic_number(case_data.atom, case_data.atomic_number)


def _test_get_atomic_number(atom, atomic_number):
    """
    Test :meth:`.Atom.get_atomic_number`

    Parameters
    ----------
    atom : :class:`.Atom`
        The atom to test.

    atomic_number : :class:`int`
        The correct atomic number of `atom`.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert atom.get_atomic_number() == atomic_number
