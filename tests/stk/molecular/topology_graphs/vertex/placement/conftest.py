import pytest
import numpy as np
import stk


@pytest.fixture(
    params=[
        stk.BuildingBlock('CCC'),
    ],
)
def building_block_0(request):
    return request.param.clone()


@pytest.fixture(
    params=[
        stk.BuildingBlock('CCBr', ['bromine']),
    ],
)
def building_block_1(request):
    return request.param.clone()


@pytest.fixture(
    params=[
        stk.BuildingBlock('BrCCCBr', ['bromine']),
    ],
)
def building_block_2(request):
    return request.param.clone()


@pytest.fixture(
    params=[
        stk.BuildingBlock('BrCC(CBr)CBr', ['bromine']),
    ],
)
def building_block_3(request):
    return request.param.clone()


@pytest.fixture(
    params=[
        stk.BuildingBlock('Brc1c(Br)cc(Br)c(Br)c1', ['bromine']),
    ],
)
def building_block_4(request):
    return request.param.clone()


@pytest.fixture(
    params=[
        stk.BuildingBlock('BrC1C(Br)C(Br)C(Br)C1Br', ['bromine']),
    ],
)
def building_block_5(request):
    return request.param.clone()


class _Graph:
    def __init__(
        self,
        placement_vertex,
        vertices,
        edges,
        position,
        alignments,
        functional_group_edges,
    ):
        self.placement_vertex = placement_vertex
        self.vertices = vertices
        self.edges = edges
        self.position = position
        self.alignments = alignments
        self.functional_group_edges = functional_group_edges


def host(position):
    vertex = stk.host_guest._HostVertexData(*position).get_vertex()
    return _Graph(
        placement_vertex=vertex,
        vertices=(vertex, ),
        edges=(),
        position=position,
        alignments={},
        functional_group_edges={},
    )


@pytest.fixture(
    params=[
        (np.array([1., 0., 0.]), np.array([1., 0., 0.])),
        (np.array([0., 1., 0.]), np.array([1., 0., 0.])),
    ]
)
def guest(request, position):
    start, stop = request.param
    vertex = stk.host_guest._GuestVertexData(
        *position,
        start=start,
        stop=stop,
    ).get_vertex()
    return _Graph(
        placement_vertex=vertex,
        vertices=(vertex, ),
        edges=(),
        position=position,
        alignments={},
        functional_group_edges={},
    )


@pytest.fixture(
    params=[
        pytest.lazy_fixture('host'),
        pytest.lazy_fixture('guest'),
        pytest.lazy_fixture('axle'),
    ],
)
def graph_0(request):
    return request.param
