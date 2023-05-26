from functools import partial

import numpy as np
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


def get_edges(vertex, angle):
    position1, position2 = get_points(vertex.get_position(), angle)
    vertex2 = stk.Vertex(1, position1)
    vertex3 = stk.Vertex(2, position2)
    yield stk.Edge(0, vertex, vertex2)
    yield stk.Edge(1, vertex, vertex3)


def get_points(position, angle):
    rotation_matrix = stk.rotation_matrix_arbitrary_axis(
        angle=angle,
        axis=np.array([0, 0, 1], dtype=np.float64),
    )
    points = (
        rotation_matrix
        @ np.array(
            [
                [-10, 0, 0],
                [10, 0, 0],
            ],
            dtype=np.float64,
        ).T
    )

    yield from (position + point for point in points.T)
