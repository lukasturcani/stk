import pytest
import numpy as np
from pytest_lazyfixture import lazy_fixture
import stk

from ...case_data import CaseData


vertices = stk.cage.vertices


@pytest.fixture(
    params=(
        lazy_fixture('cage1'),
        lazy_fixture('cage2'),
    ),
)
def case_data(request):
    return request.param


@pytest.fixture
def cage1(cls, position):
    return CaseData(
        vertex=cls(
            id=0,
            position=position,
        ),
        id=0,
        position=position,
        cell=np.array([0, 0, 0]),
    )


@pytest.fixture
def cage2(init_at_center, vertices_):
    return CaseData(
        vertex=init_at_center(
            id=0,
            vertices=vertices_,
        ),
        id=0,
        position=(
            sum(v.get_position() for v in vertices_) / len(vertices_)
        ),
        cell=np.array([0, 0, 0]),
    )


@pytest.fixture(
    params=(
        [0, 0, 0],
        [1, 2, -20],
    ),
)
def position(request):
    return np.array(request.param, dtype=np.float64)


@pytest.fixture(
    params=(
        vertices._LinearCageVertex,
        vertices._NonLinearCageVertex,
    ),
)
def cls(request):
    return request.param


@pytest.fixture(
    params=(
        vertices._LinearCageVertex.init_at_center,
        vertices._NonLinearCageVertex.init_at_center,
    ),
)
def init_at_center(request):
    return request.param


@pytest.fixture(
    params=(
        (stk.Vertex(0, [1, 2, 3]), ),
        (stk.Vertex(0, [1, 2, 3]), stk.Vertex(1, [-1, 2, -32])),
    ),
)
def vertices_(request):
    return request.param
