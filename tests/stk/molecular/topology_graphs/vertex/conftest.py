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
def graph(request, make_vertex_data):
    return request.param(make_vertex_data)


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

    def inner(x, y, z):
        return linear_vertex_data(0, 0, 0, flip=True)

    return inner


@pytest.fixture
def cage_vertex():

    def inner(x, y, z):
        return stk.cage.base._CageVertexData(
            x, y, z,
            use_bonder_placement=True,
        )

    return inner


@pytest.fixture
def one_plus_one_vertex():

    def inner(x, y, z):
        return stk.cage.three_plus_three._OnePlusOneVertexData(
            x, y, z,
            edge_normal=[1, 0, 0],
            use_bonder_placement=True,
        )

    return inner


@pytest.fixture
def cof_vertex():

    def inner(x, y, z):
        return stk.cof._COFVertexData(0, 0, 0)

    return inner


@pytest.fixture
def macrocycle_vertex():

    def inner(x, y, z):
        return stk.macrocycle._CycleVertexData(
            x, y, z,
            flip=True,
            angle=1,
        )

    return inner


@pytest.fixture
def host_vertex():

    def inner(x, y, z):
        return stk.host_guest._HostVertexData(x, y, z)

    return inner


@pytest.fixture
def guest_vertex():

    def inner(x, y, z):
        return stk.host_guest._GuestVertexData(
            x, y, z,
            start=[0, 0, 1],
            target=[1, 0, 0],
        )

    return inner


@pytest.fixture
def axle_vertex():

    def inner(x, y, z):
        return stk.rotaxane._AxleVertexData(x, y, z)

    return inner


@pytest.fixture
def cycle_vertex():

    def inner(x, y, z):
        return stk.rotaxane._CycleVertexData(x, y, z, True)

    return inner


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
def make_vertex_data(request):
    return request.param


@pytest.fixture
def vertex_data(make_vertex_data):
    return make_vertex_data(0, 0, 0)
