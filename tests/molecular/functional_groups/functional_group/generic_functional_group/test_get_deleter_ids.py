import itertools as it


def test_get_deleter_ids(generic_test_case):
    _test_get_deleter_ids(
        functional_group=generic_test_case.functional_group,
        deleters=generic_test_case.deleters,
    )


def _test_get_deleter_ids(functional_group, deleters):
    fg_atoms = it.zip_longest(
        functional_group.get_deleter_ids(),
        deleters,
    )
    for id_, atom in fg_atoms:
        assert id_ == atom.get_id()
