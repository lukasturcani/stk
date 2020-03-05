def test_get_building_block_vertices(case_data):
    assert (
        case_data.constructed_molecule.get_building_block_vertices()
        == case_data.building_block_vertices
    )
