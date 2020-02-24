import pytest
import stk
import numpy as np
from pytest_lazyfixture import lazy_fixture

from .._test_case import _TestCase


vertices = stk.molecular.topology_graphs.polymer.linear


@pytest.fixture(
    params=(
        lazy_fixture('center'),
        lazy_fixture('head'),
        lazy_fixture('tail'),
    ),
)
def test_case(request):
    return request.param


@pytest.fixture
def center(position, flip, center_building_block):
    point1 = position + [-10, 0, 0]
    point2 = position + [10, 0, 0]
    vertex = vertices._LinearVertex(0, position, flip)
    return _TestCase(
        vertex=vertex,
        edges=tuple(get_edges(vertex)),
        building_block=center_building_block,
        position=position,
        points=(point1, point2),
        nearest_points=get_nearest_points(flip, point1, point2),
        functional_group_edges=get_functional_group_edges(flip),
    )


@pytest.fixture
def head(position, flip, end_building_block):
    point1 = position + [-10, 0, 0]
    point2 = position + [10, 0, 0]

    if end_building_block.get_num_functional_groups() == 2:
        nearest_points = get_nearest_points(flip, point1, point2)
        if flip:
            functional_group_edges = {0: 1}
        else:
            functional_group_edges = {1: 1}
    else:
        nearest_points = {0: point2}
        functional_group_edges = {0: 1}

    vertex = vertices._HeadVertex(0, position, flip)
    return _TestCase(
        vertex=vertex,
        edges=(tuple(get_edges(vertex))[1], ),
        building_block=end_building_block,
        position=position,
        points=(point1, point2),
        nearest_points=nearest_points,
        functional_group_edges=functional_group_edges,
    )


@pytest.fixture
def tail(position, flip, end_building_block):
    point1 = position + [-10, 0, 0]
    point2 = position + [10, 0, 0]

    if end_building_block.get_num_functional_groups() == 2:
        nearest_points = get_nearest_points(flip, point1, point2)
        if flip:
            functional_group_edges = {1: 0}
        else:
            functional_group_edges = {0: 0}
    else:
        nearest_points = {0: point1}
        functional_group_edges = {0: 0}

    vertex = vertices._TailVertex(0, position, flip)
    return _TestCase(
        vertex=vertex,
        edges=(tuple(get_edges(vertex))[0], ),
        building_block=end_building_block,
        position=position,
        points=(point1, point2),
        nearest_points=nearest_points,
        functional_group_edges=functional_group_edges,
    )


def get_nearest_points(flip, point1, point2):
    if flip:
        return {0: point2, 1: point1}
    else:
        return {0: point1, 1: point2}


def get_functional_group_edges(flip):
    return {0: 1, 1: 0} if flip else {0: 0, 1: 1}


@pytest.fixture(
    params=(
        stk.BuildingBlock(
            smiles='BrCCNCBr',
            functional_groups=[stk.BromoFactory()],
        ),
    ),
)
def center_building_block(request):
    return request.param


@pytest.fixture(
    params=(
        stk.BuildingBlock(
            smiles='BrCC',
            functional_groups=[stk.BromoFactory()],
        ),
    ),
)
def _end_building_block(request):
    return request.param


@pytest.fixture(
    params=(
        lazy_fixture('center_building_block'),
        lazy_fixture('_end_building_block'),
    ),
)
def end_building_block(request):
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
