

def test_clone(building_block):
    clone = building_block.clone()
    fgs = it.zip_longest(clone.func_groups, building_block.func_groups)
    for fg1, fg2 in fgs:
        assert is_equivalent_fg(fg1, fg2)
