def test_get_id(case_data):
    """
    Test :meth:`.Atom.get_id`.

    Parameters
    ----------
    case_data : :class:`case_data`
        A test case. Holds the atom to test and its correct id.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert case_data.atom.get_id() == case_data.id
