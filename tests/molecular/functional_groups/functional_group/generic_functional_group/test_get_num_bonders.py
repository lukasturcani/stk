def test_get_num_bonders(generic_case_data):
    """
    Test :meth:`.GenericFunctionalGroup.get_num_bonders`.

    Parameters
    ----------
    generic_case_data : :class:`.GenericCaseData`
        The test case. Holds the functional group to test and the
        correct number of bonder atoms.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_get_num_bonders(
        functional_group=generic_case_data.functional_group,
        num_bonders=len(generic_case_data.bonders),
    )


def _test_get_num_bonders(functional_group, num_bonders):
    """
    Test :meth:`.GenericFunctionalGroup.get_num_bonders`.

    Parameters
    ----------
    functional_group : :class:`.GenericFunctionalGroup`
        The functional group to test.

    num_bonders : :class:`int`
        The correct number of bonder atoms.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert functional_group.get_num_bonders() == num_bonders
