from functools import partial

import numpy as np
import pytest
import stk
from scipy.spatial.distance import euclidean

from ....case_data import CaseData


@pytest.fixture
def linear(position, linear_aligner_edge, building_block_2):
    point1, point2 = points = (
        position + [-10, 0, 0],
        position + [10, 0, 0],
    )

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

    vertex = stk.cage.LinearVertex(
        id=0,
        position=position,
        aligner_edge=linear_aligner_edge,
    )

    return CaseData(
        vertex=vertex,
        edges=tuple(get_linear_edges(vertex)),
        building_block=building_block_2,
        position=position,
        alignment_tests={
            get_fg0_point: (point1 if linear_aligner_edge == 0 else point2),
            get_fg1_point: (point2 if linear_aligner_edge == 0 else point1),
        },
        functional_group_edges=(
            {0: 0, 1: 1} if linear_aligner_edge == 0 else {0: 1, 1: 0}
        ),
    )


def get_closest_point(points, point):
    return min(points, key=partial(euclidean, point))


def get_fg_position(id, building_block):
    functional_group = next(building_block.get_functional_groups(id))
    return building_block.get_centroid(
        atom_ids=functional_group.get_placer_ids(),
    )


def get_linear_edges(vertex):
    vertex2 = stk.Vertex(1, vertex.get_position() + [-10, 0, 0])
    vertex3 = stk.Vertex(2, vertex.get_position() + [10, 0, 0])
    yield stk.Edge(0, vertex, vertex2)
    yield stk.Edge(1, vertex, vertex3)


@pytest.fixture(params=(0, 1))
def linear_aligner_edge(request):
    return request.param


@pytest.fixture(
    params=(
        lambda: stk.BuildingBlock(
            smiles="BrCCNCBr",
            functional_groups=[stk.BromoFactory()],
        ),
    ),
)
def building_block_2(request) -> stk.BuildingBlock:
    return request.param()


@pytest.fixture(
    params=([1, 2, -20],),
)
def position(request):
    """
    The `position` of a vertex.

    """

    return np.array(request.param, dtype=np.float64)
