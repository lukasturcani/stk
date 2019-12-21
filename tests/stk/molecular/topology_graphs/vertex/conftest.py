import numpy as np
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
    params=[0., 1., 0.5, 2.],
)
def scale(request):
    return request.param


@pytest.fixture(
    params=[
        np.array([0., 0., 0.]),
        np.array([10., 3., 5.]),
        np.array([-1., -2., 0.]),
        np.array([-1., -20., -3.]),
        np.array([-23., 34., 12.]),
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


@pytest.fixture
def linear_vertex(linear_vertex_data):
    return linear_vertex_data(0, 0, 0, flip=True)


@pytest.fixture
def cage_vertex():
    return stk.cage.base._CageVertexData(
        0, 0, 0,
        use_bonder_placement=True,
    )


@pytest.fixture
def one_plus_one_vertex():
    return stk.cage.three_plus_three._OnePlusOneVertexData(
        0, 0, 0,
        edge_normal=[1, 0, 0],
        use_bonder_placement=True,
    )


@pytest.fixture
def cof_vertex():
    return stk.cof._COFVertexData(0, 0, 0)


@pytest.fixture
def macrocycle_vertex():
    return stk.macrocycle._CycleVertexData(
        0, 0, 0,
        flip=True,
        angle=1,
    )


@pytest.fixture
def host_vertex():
    return stk.host_guest._HostVertexData(0, 0, 0)


@pytest.fixture
def guest_vertex():
    return stk.host_guest._GuestVertexData(
        0, 0, 0,
        start=[0, 0, 1],
        target=[1, 0, 0],
    )


@pytest.fixture
def axle_vertex():
    return stk.rotaxane._AxleVertexData(0, 0, 0)


@pytest.fixture
def cycle_vertex():
    return stk.rotaxane._CycleVertexData(0, 0, 0, True)


@pytest.fixture(
    params=[
        pytest.lazy_fixture('linear_vertex'),
        pytest.lazy_fixture('cage_vertex'),
        pytest.lazy_fixture('one_plus_one_vertex'),
        pytest.lazy_fixture('cof_vertex'),
        pytest.lazy_fixture('macrocycle_vertex'),
        pytest.lazy_fixture('host_vertex'),
        pytest.lazy_fixture('guest_vertex'),
        pytest.lazy_fixture('axle_vertex'),
        pytest.lazy_fixture('cycle_vertex'),
    ],
)
def vertex_data(request):
    return request.param
