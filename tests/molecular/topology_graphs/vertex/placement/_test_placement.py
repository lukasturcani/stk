from scipy.spatial.distance import euclidean
import numpy as np


def test_placement(test_case):
    _test_placement(
        vertex=test_case.vertex,
        edges=test_case.edges,
        building_block=test_case.building_block,
        position=test_case.position,
        points=test_case.points,
        nearest_points=test_case.nearest_points,
        functional_group_edges=test_case.functional_group_edges,
    )


def _test_placement(
    vertex,
    edges,
    building_block,
    position,
    points,
    nearest_points,
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
    correct_points(
        points1=get_nearest_points(building_block, points),
        points2=nearest_points,
    )
    result = vertex.map_functional_groups_to_edges(
        building_block=building_block,
        edges=edges,
    )
    assert result == functional_group_edges


def get_nearest_points(building_block, points):
    distances = []
    for i, fg in enumerate(building_block.get_functional_groups()):
        fg_position = building_block.get_centroid(fg.get_placer_ids())
        for point in points:
            distances.append(
                (euclidean(fg_position, point), i, point)
            )
    nearest_points = {}
    for _, fg, point in sorted(distances, key=lambda x: x[0]):
        if fg not in nearest_points:
            nearest_points[fg] = point
    return nearest_points


def correct_points(points1, points2):
    for functional_group, point1 in points1.items():
        assert np.all(np.equal(point1, points2[functional_group]))
