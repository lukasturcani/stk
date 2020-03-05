import itertools as it


def test_get_placer_ids(case_data):
    _test_get_placer_ids(
        functional_group=case_data.functional_group,
        placers=case_data.placers,
    )


def _test_get_placer_ids(functional_group, placers):
    fg_placers = it.zip_longest(
        functional_group.get_placer_ids(),
        placers,
    )
    for id_, placer in fg_placers:
        assert id_ == placer.get_id()
