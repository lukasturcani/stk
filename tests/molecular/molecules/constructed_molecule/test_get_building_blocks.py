def test_get_building_blocks(case_data):
    assert (
        tuple(case_data.constructed_molecule.get_building_blocks())
        == case_data.building_blocks
    )
