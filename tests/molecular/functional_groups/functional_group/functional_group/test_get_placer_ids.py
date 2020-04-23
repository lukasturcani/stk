import itertools as it


def test_get_placer_ids(case_data):
    """
    Test :meth:`.FunctionalGroup.get_placer_ids`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the functional group to test and the atoms
        holding the correct ids.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_get_placer_ids(
        functional_group=case_data.functional_group,
        placers=case_data.placers,
    )


def _test_get_placer_ids(functional_group, placers):
    """
    Test :meth:`.FunctionalGroup.get_placer_ids`.

    Parameters
    ----------
    functional_group : :class:`.FunctionalGroup`
        The functional group to test.

    placers : :class:`tuple` of :class:`.Atom`
        The atoms with correct ids.

    Returns
    -------
    None : :class:`NoneType`

    """

    for id_, placer in it.zip_longest(
        functional_group.get_placer_ids(),
        placers,
    ):
        assert id_ == placer.get_id()
