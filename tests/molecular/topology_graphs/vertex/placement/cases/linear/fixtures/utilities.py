from functools import partial

import stk
from scipy.spatial.distance import euclidean


def get_fg_position(id, building_block):
    (functional_group,) = building_block.get_functional_groups(id)
    return building_block.get_centroid(
        atom_ids=functional_group.get_placer_ids(),
    )


def get_centroid(building_block):
    return building_block.get_centroid()


def get_closest_point(points, point):
    return min(points, key=partial(euclidean, point))


def get_edges(vertex):
    vertex2 = stk.Vertex(1, vertex.get_position() + [-10, 0, 0])
    vertex3 = stk.Vertex(2, vertex.get_position() + [10, 0, 0])
    yield stk.Edge(0, vertex, vertex2)
    yield stk.Edge(1, vertex, vertex3)
