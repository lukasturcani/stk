import numpy as np
from scipy.spatial.distance import euclidean


def _test_placement(graph, building_block, position):
    graph.placement_vertex.place_building_block(
        building_block=building_block,
        vertices=graph.vertices,
        edges=graph.edges,
    )
    assert (
        np.all(np.equal(building_block.get_centroid(), position))
    )


def _test_alignment(graph, building_block):
    bonder_centroids = enumerate(
        building_block.get_bonder_centroids()
    )
    for i, bonder_centroid in bonder_centroids:
        closest_edge = min(
            graph.placement_vertex.get_edge_ids(),
            key=lambda edge_id: euclidean(
                position1=bonder_centroid,
                position2=graph.edges[edge_id].get_position(),
            ),
        )
        assert closest_edge == graph.alignments[i]


def _test_assignment(graph, building_block):
    vertex = graph.placement_vertex
    functional_group_edges = vertex.assign_func_groups_to_edges(
        building_block=building_block,
        vertices=graph.vertices,
        edges=graph.edges,
    )
    assert (
        functional_group_edges == graph.functional_group_edges
    )


def _test_place_building_block(graph, building_block):
    _test_placement(graph, building_block)
    _test_alignment(graph, building_block)
    _test_assignment(graph, building_block)


def test_place_building_block_0(graph_0, building_block_0):
    _test_place_building_block(graph_0, building_block_0)


def test_place_building_block_1(graph_1, building_block_1):
    _test_place_building_block(graph_1, building_block_1)


def test_place_building_block_2(graph_2, building_block_2):
    _test_place_building_block(graph_2, building_block_2)


def test_place_building_block_3(graph_3, building_block_3):
    _test_place_building_block(graph_3, building_block_3)


def test_place_building_block_4(graph_4, building_block_4):
    _test_place_building_block(graph_4, building_block_4)


def test_place_building_block_5(graph_5, building_block_5):
    _test_place_building_block(graph_5, building_block_5)
