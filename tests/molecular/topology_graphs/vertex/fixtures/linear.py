import pytest
import stk
import numpy as np
from pytest_lazyfixture import lazy_fixture

from ._test_case import _TestCase


vertices = stk.molecular.topology_graphs.polymer.linear.vertices


@pytest.fixture(
    params=(
        lazy_fixture('center'),
        lazy_fixture('head'),
        lazy_fixture('tail'),
    ),
)
def linear(request):
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
