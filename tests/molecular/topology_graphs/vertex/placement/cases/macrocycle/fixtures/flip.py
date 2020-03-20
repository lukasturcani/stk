import pytest
from scipy.spatial.distance import euclidean
import stk
import numpy as np

from ....case_data import CaseData
from .utilities import (
    get_closest_point,
    get_fg_position,
    get_edges,
    get_points,
)


vertices = stk.molecular.topology_graphs.macrocycle.vertices


@pytest.fixture
def flip(position, building_block_2, angle):

    point2, point1 = points = tuple(get_points(position, angle))

    def get_fg0_point(building_block):
        return get_closest_point(
            points=points,
            point=get_fg_position(0, building_block),
        )

    def get_fg1_point(building_block):
        return get_closest_point(
            points=points,
            point=get_fg_position(1, building_block),
        )

    vertex = vertices._CycleVertex(0, position, True, angle+np.pi/2)
    edges = tuple(get_edges(vertex))

    fg0, = building_block_2.get_functional_groups(0)
    fg0_position = building_block_2.get_centroid(fg0.get_placer_ids())

    def fg0_distance(edge):
        return euclidean(fg0_position, edge.get_position())

    edge0 = min(edges, key=fg0_distance)
    edge1 = edges[1] if edge0 is edges[0] else edges[0]

    return CaseData(
        vertex=vertex,
        edges=edges,
        building_block=building_block_2,
        position=position,
        alignment_tests={
            get_fg0_point: point1,
            get_fg1_point: point2,
        },
        functional_group_edges={0: edge0.get_id(), 1: edge1.get_id()},
    )
