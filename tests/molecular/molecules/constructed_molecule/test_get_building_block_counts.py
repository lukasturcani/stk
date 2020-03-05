def test_get_building_block_counts(case_data):
    assert (
        case_data.constructed_molecule.get_building_block_counts()
        == case_data.building_block_counts
    )
