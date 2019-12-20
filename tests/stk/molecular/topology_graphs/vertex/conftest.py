import numpy as np
import pytest
import fixtures
import stk


@pytest.fixture(
    params=[
        fn for name, fn in fixtures.make_vertex.__dict__.items()
        if not name.startswith('_') and callable(fn)
    ],
)
def make_vertex(request):
    return request.param


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
    params=[
        (0, 3, 5),
        (5, 3, 2),
        (1, 12),
        (),
    ],
)
def edge_ids(request):
    return tuple(request.param)


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
