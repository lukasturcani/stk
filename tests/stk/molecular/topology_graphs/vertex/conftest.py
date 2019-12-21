import numpy as np
import pytest_cases
import pytest
import fixtures
import stk


@pytest.fixture(
    params=[
        fn for name, fn in fixtures.graph.__dict__.items()
        if not name.startswith('_') and callable(fn)
    ],
)
def graph(request):
    return request.param()


@pytest.fixture(
    params=[
        lambda molecule: None,
        lambda molecule: range(0, len(molecule.atoms), 2),
    ],
)
def get_atom_ids(request):
    return request.param


@pytest.fixture(
    params=[
        stk.ConstructedMolecule(
            building_blocks=[stk.BuildingBlock('BrCCBr', ['bromine'])],
            topology_graph=stk.polymer.Linear('A', 3),
        ),
    ]
)
def constructed_molecule(request):
    return request.param.clone()


@pytest.fixture(
    params=[0, 1, 0.5, 2],
)
def scale(request):
    return request.param


@pytest.fixture(
    params=[
        np.array([0, 0, 0]),
        np.array([10, 3, 5]),
        np.array([-1, -2, 0]),
        np.array([-1, -20, -3]),
        np.array([-23, 34, 12]),
    ],
)
def position(request):
    return np.array(request.param)


@pytest.fixture(
    params=[0, 1, 2, 3],
)
def add_edge_data(num_edges):

    def inner(data):
        edges = [
            stk.EdgeData(data, position=[1, 2, 3])
            for i in range(num_edges)
        ]
        for i, edge in enumerate(edges):
            edge.id = i
        return edges

    return inner


@pytest.fixture(
    params=[
        (0, 0, 0),
        (1, 0, 0),
        (-1, 0, 0),
        (10, 3, -1),
    ],
)
def cell(request):
    return tuple(request.param)


@pytest.fixture(
    params=[
        stk.polymer._LinearVertexData,
        stk.polymer._TailVertexData,
        stk.polymer._HeadVertexData,
    ]
)
def linear_vertex_data(request):
    return request.param


@pytest_cases.pytest_fixture_plus
def linear_vertex(
    linear_vertex_data,
    position,
    add_edge_data,
    cell,
    flip
):
    data = linear_vertex_data(
        position=position,
        cell=cell,
        flip=flip,
    )
    add_edge_data(data)
    return data
