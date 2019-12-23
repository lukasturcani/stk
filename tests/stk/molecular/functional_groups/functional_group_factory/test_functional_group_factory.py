import itertools as it


def _test_equivalent_atoms(atoms1, atoms2):
    for a1, a2 in it.zip_longest(atoms1, atoms2):
        assert a1 is not a2
        assert a1.id == a2.id
        assert a1.__class__ is a2.__class__


def _test_equivalent_functional_group(
    functional_group1,
    functional_group2,
):
    _test_equivalent_atoms(
        atoms1=functional_group1.get_atoms(),
        atoms2=functional_group2.get_atoms(),
    )
    _test_equivalent_atoms(
        atoms1=functional_group1.get_bonders(),
        atoms2=functional_group2.get_bonders(),
    )
    _test_equivalent_atoms(
        atoms1=functional_group1.get_deleters(),
        atoms2=functional_group2.get_deleters(),
    )


def test_get_functional_groups(get_functional_groups_test_case):
    fgs = it.zip_longest(
        get_functional_groups_test_case.functional_groups,
        get_functional_groups_test_case.factory.get_functional_groups(
            molecule=get_functional_groups_test_case.molecule,
        ),
    )
    for expected_fg, fg in fgs:
        _test_equivalent_functional_group(expected_fg, fg)
