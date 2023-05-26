from collections import Counter


def check_building_block_counts(old_state, new_state, building_blocks):
    counts = Counter(building_blocks)
    old_count = old_state.get_building_block_counts()
    new_count = new_state.get_building_block_counts()

    for building_block in old_count:
        expected_count = old_count[building_block] - new_count[building_block]
        assert counts.get(building_block, 0) == expected_count

    new_only = new_count.keys() - old_count.keys()
    added_only = counts.keys() - old_count.keys()
    assert len(new_only) == len(added_only)

    for building_block in new_only:
        assert new_count[building_block] == counts[building_block]
