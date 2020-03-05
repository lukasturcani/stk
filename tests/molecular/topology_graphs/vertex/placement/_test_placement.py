import numpy as np


def test_placement(case_data):
    _test_placement(
        vertex=case_data.vertex,
        edges=case_data.edges,
        building_block=case_data.building_block,
        position=case_data.position,
        alignment_tests=case_data.alignment_tests,
        functional_group_edges=case_data.functional_group_edges,
    )


def _test_placement(
    vertex,
    edges,
    building_block,
    position,
    alignment_tests,
    functional_group_edges,
):
    position_matrix = vertex.place_building_block(
        building_block=building_block,
        edges=edges,
    )
    building_block = building_block.with_position_matrix(
        position_matrix=position_matrix,
    )
    assert np.allclose(
        a=building_block.get_centroid(building_block.get_placer_ids()),
        b=position,
        atol=1e-14,
    )
    for test, result in alignment_tests.items():
        assert np.all(np.equal(test(building_block), result))

    functional_group_edges_ = vertex.map_functional_groups_to_edges(
        building_block=building_block,
        edges=edges,
    )
    assert functional_group_edges_ == functional_group_edges
