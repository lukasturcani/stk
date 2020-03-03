import stk
import numpy as np
import pytest

from ._test_case import _TestCase


vertices = stk.molecular.topology_graphs.cof.vertices


@pytest.fixture
def cof(cls, id, position, aligner_edge, cell):
    return _TestCase(
        vertex=cls(
            id=id,
            position=position,
            aligner_edge=aligner_edge,
            cell=cell,
        ),
        id=id,
        cell=cell,
        position=position,
    )


@pytest.fixture(
    params=(
        vertices._LinearCofVertex,
        vertices._NonLinearCofVertex,
    ),
)
def cls(request):
    return request.param


@pytest.fixture(params=(0, 20))
def id(request):
    return request.param


@pytest.fixture(params=(0, 1))
def aligner_edge(request):
    return request.param


@pytest.fixture(
    params=(
        np.array([0, 0, 0]),
        np.array([-20, 1, 21]),
    ),
)
def cell(request):
    return np.array(request.param)
