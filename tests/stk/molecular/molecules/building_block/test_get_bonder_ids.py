

def test_get_bonder_ids(building_block, get_fg_ids):
    fg_ids = get_fg_ids(building_block)
    if fg_ids is None:
        fg_ids = range(building_block.get_num_functional_groups)

    expected_bonder_ids = (
        bid
        for fg in building_block.get_functional_groups(fg_ids)
        for bid in fg.get_bonder_ids()
    )
    ids = building_block.get_bonder_ids(
        fg_ids=get_fg_ids(building_block),
    )
    bonder_ids = it.zip_longest(ids, expected_bonder_ids)
    for bonder_id, expected_bonder_id in bonder_ids:
        assert bonder_id == expected_bonder_id
