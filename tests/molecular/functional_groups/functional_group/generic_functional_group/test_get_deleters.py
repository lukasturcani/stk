import itertools as it


def test_get_deleters(generic_case_data):
    """
    Test :meth:`.GenericFunctionalGroup.get_deleters`.

    Parameters
    ----------
    generic_case_data : :class:`.GenericCaseData`
        The test case. Holds the functional group to test and the
        correct deleter atoms.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_get_deleters(
        functional_group=generic_case_data.functional_group,
        deleters=generic_case_data.deleters,
    )


def _test_get_deleters(functional_group, deleters):
    """
    Test :meth:`.GenericFunctionalGroup.get_deleters`.

    Parameters
    ----------
    functional_group : :class:`.GenericFunctionalGroup`
        The functional group to test.

    deleters : :class:`tuple` of :class:`.Atom`
        The correct deleter atoms.

    Returns
    -------
    None : :class:`NoneType`

    """

    for atom1, atom2 in it.zip_longest(
        functional_group.get_deleters(),
        deleters,
    ):
        assert atom1 is atom2
