import pytest
import stk
import numpy as np
from pytest_lazyfixture import lazy_fixture
from scipy.spatial.distance import euclidean
from functools import partial

from .._test_case import _TestCase


vertices = stk.molecular.topology_graphs.polymer.linear


@pytest.fixture(
    params=(
        lazy_fixture('center'),
        lazy_fixture('head_1'),
        lazy_fixture('head_2'),
        lazy_fixture('head_3'),
        lazy_fixture('tail_1'),
        lazy_fixture('tail_2'),
        lazy_fixture('tail_3'),
    ),
)
def test_case(request):
    return request.param


def get_fg_position(id, building_block):
    functional_group = next(building_block.get_functional_groups(id))
    return building_block.get_centroid(
        atom_ids=functional_group.get_placer_ids(),
    )


def get_centroid(building_block):
    return building_block.get_centroid()


def get_closest_point(points, point):
    return min(points, key=partial(euclidean, point))


@pytest.fixture
def center(position, flip, building_block_2):

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

    vertex = vertices._LinearVertex(0, position, flip)
    return _TestCase(
        vertex=vertex,
        edges=tuple(get_edges(vertex)),
        building_block=building_block_2,
        position=position,
        alignment_tests={
            get_fg0_point: point2 if flip else point1,
            get_fg1_point: point1 if flip else point2,
        },
        functional_group_edges={0: 1, 1: 0} if flip else {0: 0, 1: 1},
    )


@pytest.fixture
def head_1(position, flip, building_block_1):
    point1, point2 = points = (
        position + [-10, 0, 0],
        position + [10, 0, 0],
    )

    def get_centroid_point(building_block):
        return get_closest_point(
            points=points,
            point=get_centroid(building_block),
        )

    vertex = vertices._HeadVertex(0, position, flip)
    return _TestCase(
        vertex=vertex,
        edges=(tuple(get_edges(vertex))[1], ),
        building_block=building_block_1,
        position=position,
        alignment_tests={get_centroid_point: point1},
        functional_group_edges={0: 1},
    )


@pytest.fixture
def head_2(position, building_block_2):
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

    vertex = vertices._HeadVertex(0, position, False)
    return _TestCase(
        vertex=vertex,
        edges=(tuple(get_edges(vertex))[1], ),
        building_block=building_block_2,
        position=position,
        alignment_tests={
            get_fg0_point: point1,
            get_fg1_point: point2,
        },
        functional_group_edges={1: 1},
    )


@pytest.fixture
def head_3(position, building_block_2):
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

    vertex = vertices._HeadVertex(0, position, True)
    return _TestCase(
        vertex=vertex,
        edges=(tuple(get_edges(vertex))[1], ),
        building_block=building_block_2,
        position=position,
        alignment_tests={
            get_fg0_point: point2,
            get_fg1_point: point1,
        },
        functional_group_edges={0: 1},
    )


@pytest.fixture
def tail_1(position, flip, building_block_1):
    point1, point2 = points = (
        position + [-10, 0, 0],
        position + [10, 0, 0],
    )

    def get_centroid_point(building_block):
        return get_closest_point(
            points=points,
            point=get_centroid(building_block),
        )

    vertex = vertices._TailVertex(0, position, flip)
    return _TestCase(
        vertex=vertex,
        edges=(tuple(get_edges(vertex))[0], ),
        building_block=building_block_1,
        position=position,
        alignment_tests={get_centroid_point: point2},
        functional_group_edges={0: 0},
    )


@pytest.fixture
def tail_2(position, building_block_2):
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

    vertex = vertices._TailVertex(0, position, False)
    return _TestCase(
        vertex=vertex,
        edges=(tuple(get_edges(vertex))[0], ),
        building_block=building_block_2,
        position=position,
        alignment_tests={
            get_fg0_point: point1,
            get_fg1_point: point2,
        },
        functional_group_edges={0: 0},
    )


@pytest.fixture
def tail_3(position, building_block_2):
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

    vertex = vertices._TailVertex(0, position, True)
    return _TestCase(
        vertex=vertex,
        edges=(tuple(get_edges(vertex))[0], ),
        building_block=building_block_2,
        position=position,
        alignment_tests={
            get_fg0_point: point2,
            get_fg1_point: point1,
        },
        functional_group_edges={1: 0},
    )


@pytest.fixture(
    params=(
        stk.BuildingBlock(
            smiles='BrCCNCBr',
            functional_groups=[stk.BromoFactory()],
        ),
    ),
)
def building_block_2(request):
    return request.param


@pytest.fixture(
    params=(
        stk.BuildingBlock(
            smiles='BrCC',
            functional_groups=[stk.BromoFactory()],
        ),
    ),
)
def building_block_1(request):
    return request.param


@pytest.fixture(params=(True, False))
def flip(request):
    return request.param


@pytest.fixture(
    params=(
        [0, 0, 0],
        [1, 2, -20],
    ),
)
def position(request):
    return np.array(request.param, dtype=np.float64)


def get_edges(vertex):
    vertex2 = stk.Vertex(1, vertex.get_position() + [-10, 0, 0])
    vertex3 = stk.Vertex(2, vertex.get_position() + [10, 0, 0])
    yield stk.Edge(0, vertex, vertex2)
    yield stk.Edge(1, vertex, vertex3)
