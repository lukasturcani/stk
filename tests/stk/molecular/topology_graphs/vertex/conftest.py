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
def add_edge_data(request):

    def inner(data):
        edges = [
            stk.EdgeData(data, position=[1, 2, 3])
            for i in range(request.param)
        ]
        for i, edge in enumerate(edges):
            edge.id = i
        return data

    return inner


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
):
    data = linear_vertex_data(
        position=position,
        flip=True,
    )
    return add_edge_data(data)


@pytest_cases.pytest_fixture_plus
def cage_vertex(position, add_edge_data):
    data = stk.cage._CageVertexData(
        *position,
        use_bonder_placement=True,
    )
    return add_edge_data(data)


@pytest_cases.pytest_fixture_plus
def one_plus_one_vertex(position, add_edge_data):
    data = stk.cage._OnePlusOneVertexData(
        *position,
        edge_normal=[1, 0, 0],
        use_bonder_placement=True,
    )
    return add_edge_data(data)


@pytest_cases.pytest_fixture_plus
def cof_vertex(position, add_edge_data):
    data = stk.cof._COFVertexData(*position)
    return add_edge_data(data)


@pytest_cases.pytest_fixture_plus
def macrocycle_vertex(position, add_edge_data):
    data = stk.macrocycle._CycleVertexData(
        *position,
        flip=True,
        angle=1,
    )
    return add_edge_data(data)


@pytest_cases.pytest_fixture_plus
def host_vertex(position, add_edge_data):
    data = stk.host_guest._HostVertexData(*position)
    return add_edge_data(data)


@pytest_cases.pytest_fixture_plus
def guest_vertex(position, add_edge_data):
    data = stk.host_guest._GuestVertexData(
        *position,
        start=[0, 0, 1],
        target=[1, 0, 0],
    )
    return add_edge_data(data)


@pytest_cases.pytest_fixture_plus
def axle_vertex(position, add_edge_data):
    data = stk.rotaxane._AxleVertexData(*position)
    return add_edge_data(data)


@pytest_cases.pytest_fixture_plus
def cycle_vertex(position, add_edge_data):
    data = stk.rotaxane._CycleVertexData(*position, True)
    return add_edge_data(data)


pytest_cases.fixture_union(
    name='vertex_data',
    fixtures=[
        linear_vertex,
        cage_vertex,
        one_plus_one_vertex,
        cof_vertex,
        macrocycle_vertex,
        host_vertex,
        guest_vertex,
        axle_vertex,
        cycle_vertex,
    ]
)
