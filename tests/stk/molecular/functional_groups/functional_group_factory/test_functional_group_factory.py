import itertools as it


def _test_equivalent_atoms(atoms1, atoms2):
    for a1, a2 in it.zip_longest(atoms1, atoms2):
        assert a1 is not a2
        assert a1.id == a2.id
        assert a1.__class__ is a2.__class__


def atom_id(atom):
    return atom.id


def _test_equivalent_functional_group(
    functional_group1,
    functional_group2,
):
    assert functional_group1.__class__ is functional_group2.__class__

    _test_equivalent_atoms(
        atoms1=sorted(functional_group1.get_atoms(), key=atom_id),
        atoms2=sorted(functional_group2.get_atoms(), key=atom_id),
    )
    _test_equivalent_atoms(
        atoms1=sorted(functional_group1.get_bonders(), key=atom_id),
        atoms2=sorted(functional_group2.get_bonders(), key=atom_id),
    )
    _test_equivalent_atoms(
        atoms1=sorted(functional_group1.get_deleters(), key=atom_id),
        atoms2=sorted(functional_group2.get_deleters(), key=atom_id),
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
