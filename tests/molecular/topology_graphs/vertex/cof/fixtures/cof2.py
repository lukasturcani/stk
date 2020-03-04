import pytest
import stk

from ..._test_case import _TestCase

vertices = stk.molecular.topology_graphs.cof.vertices


@pytest.fixture
def cof2(init_at_center, id, vertices, aligner_edge, cell):
    return _TestCase(
        vertex=init_at_center(id, vertices, aligner_edge, cell),
        id=id,
        position=sum(v.get_position() for v in vertices)/len(vertices),
        cell=cell,
    )


@pytest.fixture(
    params=(
        vertices._LinearCofVertex.init_at_center,
        vertices._NonLinearCofVertex.init_at_center,
    ),
)
def init_at_center(request):
    return request.param
