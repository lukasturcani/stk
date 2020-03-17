import itertools as it


def test_get_atom_ids(case_data):
    """
    Test :meth:`.FunctionalGroup.get_atom_ids`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        The test case. Holds the functional group to test and the
        correct atom ids.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_get_atom_ids(
        functional_group=case_data.functional_group,
        atoms=case_data.atoms,
    )


def _test_get_atom_ids(functional_group, atoms):
    """
    Test :meth:`.FunctionalGroup.get_atom_ids`.

    Parameters
    ----------
    functional_group : :class:`.FunctionalGroup`
        The functional group to test.

    atoms : :class:`tuple` of :class:`.Atom`
        The atoms with correct ids.

    Returns
    -------
    None : :class:`NoneType`

    """

    for id_, atom in it.zip_longest(
        functional_group.get_atom_ids(),
        atoms,
    ):
        assert id_ == atom.get_id()
