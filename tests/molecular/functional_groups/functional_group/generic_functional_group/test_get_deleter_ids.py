import itertools as it


def test_get_deleter_ids(generic_case_data):
    """
    Test :meth:`.GenericFunctionalGroup.get_deleter_ids`.

    Parameters
    ----------
    generic_case_data : :class:`.GenericCaseData`
        The test case. Holds the functional group to test and the
        atoms holding the correct deleter ids.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_get_deleter_ids(
        functional_group=generic_case_data.functional_group,
        deleters=generic_case_data.deleters,
    )


def _test_get_deleter_ids(functional_group, deleters):
    """
    Test :meth:`.GenericFunctionalGroup.get_deleter_ids`.

    Parameters
    ----------
    functional_group : :class:`.GenericFunctionalGroup`
        The functional group to test.

    deleters : :class:`tuple` of :class:`.Atom`
        The atoms holding the correct deleter ids.

    Returns
    -------
    None : :class:`NoneType`

    """

    for id_, atom in it.zip_longest(
        functional_group.get_deleter_ids(),
        deleters,
    ):
        assert id_ == atom.get_id()
