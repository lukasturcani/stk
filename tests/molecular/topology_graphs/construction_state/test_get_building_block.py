import itertools as it


def test_get_building_block(test_case):
    _test_get_building_block(
        construction_state=test_case.construction_state,
        building_blocks=test_case.building_blocks,
    )


def _test_get_building_block(construction_state, building_blocks):
    state_bbs = map(
        construction_state.get_building_block,
        range(construction_state.get_num_vertices()),
    )
    building_blocks_ = it.zip_longest(state_bbs, building_blocks)
    for bb1, bb2 in building_blocks_:
        assert bb1 is bb2
