import itertools as it


def test_get_placer_ids(test_case):
    _test_get_placer_ids(
        functional_group=test_case.functional_group,
        placers=test_case.placers,
    )


def _test_get_placer_ids(functional_group, placers):
    fg_placers = it.zip_longest(
        functional_group.get_placer_ids(),
        placers,
    )
    for id_, placer in fg_placers:
        assert id_ == placer.get_id()
