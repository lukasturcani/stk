import itertools as it


def test_get_bonders(generic_test_case):
    _test_get_bonders(
        functional_group=generic_test_case.functional_group,
        bonders=generic_test_case.bonders,
    )


def _test_get_bonders(functional_group, bonders):
    fg_atoms = it.zip_longest(functional_group.get_bonders(), bonders)
    for atom1, atom2 in fg_atoms:
        atom1 is atom2
