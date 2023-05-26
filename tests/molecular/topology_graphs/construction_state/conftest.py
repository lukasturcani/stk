import pytest
import stk


@pytest.fixture
def construction_state(building_block_vertices):
    return stk.ConstructionState(
        building_block_vertices=building_block_vertices,
        edges=get_edges(building_block_vertices),
        scale=1,
    )


def get_edges(building_block_vertices):
    vertices = (
        vertex
        for vertices in building_block_vertices.values()
        for vertex in vertices
    )
    vertex1 = next(vertices)
    for id_, vertex2 in enumerate(vertices):
        yield stk.Edge(id_, vertex1, vertex2)
        vertex1 = vertex2


@pytest.fixture
def building_block_vertices(building_block1, building_block2):
    return {
        building_block1: (
            stk.Vertex(0, [0, 0, 0]),
            stk.Vertex(1, [10, 0, 0]),
        ),
        building_block2: (
            stk.Vertex(2, [20, 0, 0]),
            stk.Vertex(3, [30, 0, 0]),
        ),
    }


@pytest.fixture(
    params=(stk.BuildingBlock("BrCC", [stk.BromoFactory()]),),
)
def building_block1(request):
    return request.param


@pytest.fixture(
    params=(stk.BuildingBlock("BrCCBr", [stk.BromoFactory()]),),
)
def building_block2(request):
    return request.param
