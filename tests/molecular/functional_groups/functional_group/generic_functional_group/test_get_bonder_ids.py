import itertools as it


def test_get_bonder_ids(generic_case_data):
    """
    Test :meth:`.GenericFunctionalGroup.get_bonder_ids`.

    Parameters
    ----------
    generic_case_data : :class:`.GenericCaseData`
        The test case. Holds the functional group to test and the
        atoms holding the correct bonder ids.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_get_bonder_ids(
        functional_group=generic_case_data.functional_group,
        bonders=generic_case_data.bonders,
    )


def _test_get_bonder_ids(functional_group, bonders):
    """
    Test :meth:`.GenericFunctionalGroup.get_bonder_ids`.

    Parameters
    ----------
    functional_group : :class:`.GenericFunctionalGroup`
        The functional group to test.

    bonders : :class:`tuple` of :class:`.Atom`
        The atoms holding the correct bonder ids.

    Returns
    -------
    None : :class:`NoneType`

    """

    for id_, atom in it.zip_longest(
        functional_group.get_bonder_ids(),
        bonders,
    ):
        assert id_ == atom.get_id()
