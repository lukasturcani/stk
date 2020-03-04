import pytest
import stk

from ..._test_case import _TestCase

vertices = stk.molecular.topology_graphs.cof.vertices


@pytest.fixture
def cof1(cls, id, position, aligner_edge, cell):
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
