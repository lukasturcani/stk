import itertools as it


def test_get_deleter_ids(generic_case_data):
    _test_get_deleter_ids(
        functional_group=generic_case_data.functional_group,
        deleters=generic_case_data.deleters,
    )


def _test_get_deleter_ids(functional_group, deleters):
    fg_atoms = it.zip_longest(
        functional_group.get_deleter_ids(),
        deleters,
    )
    for id_, atom in fg_atoms:
        assert id_ == atom.get_id()
