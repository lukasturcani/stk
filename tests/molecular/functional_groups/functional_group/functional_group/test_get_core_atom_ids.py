import itertools as it


def test_get_core_atom_ids(case_data):
    """
    Test :meth:`.FunctionalGroup.get_core_atom_ids`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        The test case. Holds the functional group to test and the
        atoms holding the correct ids.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_get_core_atom_ids(
        functional_group=case_data.functional_group,
        core_atoms=case_data.core_atoms,
    )


def _test_get_core_atom_ids(functional_group, core_atoms):
    """
    Test :meth:`.FunctionalGroup.get_core_atom_ids`.

    Parameters
    ----------
    functional_group : :class:`.FunctionalGroup`
        The functional group to test.

    core_atoms : :class:`tuple` of :class:`.Atom`
        Atoms with the correct ids.

    Returns
    -------
    None : :class:`NoneType`

    """

    for id_, core_atom in it.zip_longest(
        functional_group.get_core_atom_ids(),
        core_atoms,
    ):
        assert id_ == core_atom.get_id()
