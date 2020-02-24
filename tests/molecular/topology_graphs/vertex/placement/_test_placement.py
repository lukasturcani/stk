from scipy.spatial.distance import euclidean
import numpy as np


def test_placement(test_case):
    _test_placement(
        vertex=test_case.vertex,
        edges=test_case.edges,
        building_block=test_case.building_block,
        position=test_case.position,
        points=test_case.points,
        functional_group_edges=test_case.functional_group_edges,
    )


def _test_placement(
    vertex,
    edges,
    building_block,
    position,
    points,
    functional_group_edges,
):
    position_matrix = vertex.place_building_block(building_block)
    building_block = building_block.with_position_matrix(
        position_matrix=position_matrix,
    )
    assert np.allclose(
        a=building_block.get_centroid(building_block.get_placer_ids()),
        b=position,
        atol=1e-14,
    )
    for get_point, point in points.items():
        assert np.all(np.equal(get_point(building_block), point))

    functional_group_edges_ = vertex.map_functional_groups_to_edges(
        building_block=building_block,
        edges=edges,
    )
    assert functional_group_edges_ == functional_group_edges
