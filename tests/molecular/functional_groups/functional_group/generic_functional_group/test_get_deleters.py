import itertools as it


def test_get_deleters(generic_case_data):
    _test_get_deleters(
        functional_group=generic_case_data.functional_group,
        deleters=generic_case_data.deleters,
    )


def _test_get_deleters(functional_group, deleters):
    fg_atoms = it.zip_longest(
        functional_group.get_deleters(),
        deleters,
    )
    for atom1, atom2 in fg_atoms:
        atom1 is atom2
