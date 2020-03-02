def test_get_building_block_counts(test_case):
    _test_get_building_block_counts(
        construction_state=test_case.construction_state,
        building_block_counts=test_case.building_block_counts,
    )


def _test_get_building_block_counts(
    construction_state,
    building_block_counts,
):

    counts1 = construction_state.get_building_block_counts()
    all_keys = counts1.keys() | building_block_counts.keys()
    assert counts1.keys() == all_keys

    for building_block, count in counts1.items():
        assert building_block_counts[building_block] == count
