import itertools as it


def test_get_atoms(case_data):
    """
    Test :meth:`.FunctionalGroup.get_atoms`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        The test case. Holds the functional group to test and the
        correct atoms.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_get_atoms(
        functional_group=case_data.functional_group,
        atoms=case_data.atoms,
    )


def _test_get_atoms(functional_group, atoms):
    """
    Test :meth:`.FunctionalGroup.get_atoms`.

    Parameters
    ----------
    functional_group : :class:`.FunctionalGroup`
        The functional group to test.

    atoms : :class:`tuple` of :class:`.Atom`
        The correct atoms.

    Returns
    -------
    None : :class:`NoneType`

    """

    for atom1, atom2 in it.zip_longest(functional_group.get_atoms(), atoms):
        assert atom1 is atom2
