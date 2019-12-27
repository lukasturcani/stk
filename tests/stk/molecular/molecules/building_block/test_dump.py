

def test_dump_and_load(tmpdir, building_block):
    path = str(tmpdir / 'building_block.json')
    building_block.dump(path)
    new_building_block = stk.BuildingBlock.load(path)
    assert is_equivalent_building_block(
        building_block1=building_block,
        building_block2=new_building_block,
    )
