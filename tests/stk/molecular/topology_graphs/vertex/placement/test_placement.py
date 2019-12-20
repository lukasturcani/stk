import numpy as np
from scipy.spatial.distance import euclidean


def _test_placement(test_case, building_block, position):
    test_case.vertex.place_building_block(
        building_block=building_block,
        vertices=test_case.vertices,
        edges=test_case.edges,
    )
    assert (
        np.all(np.equal(building_block.get_centroid(), position))
    )


def _test_alignment(test_case, building_block):
    bonder_centroids = enumerate(
        building_block.get_bonder_centroids()
    )
    for i, bonder_centroid in bonder_centroids:
        closest_edge = min(
            test_case.vertex.get_edge_ids(),
            key=lambda edge_id: euclidean(
                position1=bonder_centroid,
                position2=test_case.edges[edge_id].get_position(),
            ),
        )
        assert closest_edge == test_case.alignments[i]


def _test_assignment(test_case, building_block):
    vertex = test_case.vertex
    functional_group_edges = vertex.assign_func_groups_to_edges(
        building_block=building_block,
        vertices=test_case.vertices,
        edges=test_case.edges,
    )
    assert (
        functional_group_edges == test_case.functional_group_edges
    )


def _test_place_building_block(
    make_test_case,
    building_block,
    position,
):
    test_case = make_test_case(position)
    _test_placement(test_case, building_block, position)
    _test_alignment(test_case, building_block)
    _test_assignment(test_case, building_block)


def test_place_building_block_0(
    make_test_case_0,
    building_block_0,
    position,
):
    _test_place_building_block(
        make_test_case=make_test_case_0,
        building_block=building_block_0,
        position=position,
    )


def test_place_building_block_1(
    make_test_case_1,
    building_block_1,
    position,
):
    _test_place_building_block(
        make_test_case=make_test_case_1,
        building_block=building_block_1,
        position=position,
    )


def test_place_building_block_2(
    make_test_case_2,
    building_block_2,
    position,
):
    _test_place_building_block(
        make_test_case=make_test_case_2,
        building_block=building_block_2,
        position=position,
    )


def test_place_building_block_3(
    make_test_case_3,
    building_block_3,
    position,
):
    _test_place_building_block(
        make_test_case=make_test_case_3,
        building_block=building_block_3,
        position=position,
    )


def test_place_building_block_4(
    make_test_case_4,
    building_block_4,
    position,
):
    _test_place_building_block(
        make_test_case=make_test_case_4,
        building_block=building_block_4,
        position=position,
    )


def test_place_building_block_5(
    make_test_case_5,
    building_block_5,
    position,
):

    _test_place_building_block(
        make_test_case=make_test_case_5,
        building_block=building_block_5,
        position=position,
    )
