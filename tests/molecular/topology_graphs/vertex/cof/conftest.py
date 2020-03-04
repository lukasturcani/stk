import stk
import numpy as np
import pytest

from ._test_case import _TestCase


vertices = stk.molecular.topology_graphs.cof.vertices


@pytest.fixture
def cof(cof_cls, id, position, aligner_edge, cell):
    return _TestCase(
        vertex=cof_cls(
            id=id,
            position=position,
            aligner_edge=aligner_edge,
            cell=cell,
        ),
        id=id,
        cell=cell,
        position=position,
    )


@pytest.fixture
def cof_center(cof_init_at_center, id, vertices, aligner_edge, cell):
    return _TestCase(
        vertex=cof_init_at_center(id, vertices, aligner_edge, cell),
        id=id,
        position=sum(v.get_position() for v in vertices)/len(vertices),
        aligner_edge=aligner_edge,
        cell=cell,
    )


@pytest.fixture
def


@pytest.fixture(
    params=(
        vertices._LinearCofVertex,
        vertices._NonLinearCofVertex,
    ),
)
def cof_cls(request):
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
