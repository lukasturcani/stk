import itertools as it


def test_get_bonder_ids(test_case):
    _test_get_bonder_ids(
        functional_group=test_case.functional_group,
        bonders=test_case.bonders,
    )


def _test_get_bonder_ids(functional_group, bonders):
    fg_atoms = it.zip_longest(
        functional_group.get_bonder_ids(),
        bonders,
    )
    for id_, atom in fg_atoms:
        assert id_ == atom.id
