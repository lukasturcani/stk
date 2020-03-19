import itertools as it


def test_get_bonders(generic_case_data):
    """
    Test :meth:`.GenericFunctionalGroup.get_bonders`.

    Parameters
    ----------
    generic_case_data : :class:`.GenericCaseData`
        The test case. Holds the functional group to test and the
        correct bonder atoms.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_get_bonders(
        functional_group=generic_case_data.functional_group,
        bonders=generic_case_data.bonders,
    )


def _test_get_bonders(functional_group, bonders):
    """
    Test :meth:`.GenericFunctionalGroup.get_bonders`.

    Parameters
    ----------
    functional_group : :class:`.GenericFunctionalGroup`
        The functional group to test.

    bonders : :class:`tuple` of :class:`.Atom`
        The correct bonder atoms.

    Returns
    -------
    None : :class:`NoneType`

    """

    for atom1, atom2 in it.zip_longest(
        functional_group.get_bonders(),
        bonders,
    ):
        assert atom1 is atom2
