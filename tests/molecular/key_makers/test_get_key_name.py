def test_get_key_name(case_data):
    """
    Test :meth:`.get_key_name`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the key maker to test and the correct name.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert case_data.key_maker.get_key_name() == case_data.key_name
