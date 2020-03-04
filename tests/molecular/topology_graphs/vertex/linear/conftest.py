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
def center(id, position, flip):
    return _TestCase(
        vertex=vertices._LinearVertex(id, position, flip),
        id=id,
        position=position,
        cell=np.array([0, 0, 0]),
    )


@pytest.fixture
def head(id, position, flip):
    return _TestCase(
        vertex=vertices._HeadVertex(id, position, flip),
        id=id,
        position=position,
        cell=np.array([0, 0, 0]),
    )


@pytest.fixture
def tail(id, position, flip):
    return _TestCase(
        vertex=vertices._TailVertex(id, position, flip),
        id=id,
        position=position,
        cell=np.array([0, 0, 0]),
    )


@pytest.fixture(params=(0, ))
def id(request):
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
