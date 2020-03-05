import itertools as it


def test_get_bonder_ids(generic_case_data):
    _test_get_bonder_ids(
        functional_group=generic_case_data.functional_group,
        bonders=generic_case_data.bonders,
    )


def _test_get_bonder_ids(functional_group, bonders):
    fg_atoms = it.zip_longest(
        functional_group.get_bonder_ids(),
        bonders,
    )
    for id_, atom in fg_atoms:
        assert id_ == atom.get_id()
