import pytest
import stk

from ....case_data import CaseData


@pytest.fixture
def cof2(init_at_center, id, vertices_, aligner_edge, cell):
    return CaseData(
        vertex=init_at_center(id, vertices_, aligner_edge, cell),
        id=id,
        position=(sum(v.get_position() for v in vertices_) / len(vertices_)),
        cell=cell,
    )


@pytest.fixture(
    params=(
        stk.cof.LinearVertex.init_at_center,
        stk.cof.NonLinearVertex.init_at_center,
        stk.cof.UnaligningVertex.init_at_center,
    ),
)
def init_at_center(request):
    return request.param


@pytest.fixture(
    params=(
        lambda: (stk.Vertex(0, [1, 2, 3]),),
        lambda: (
            stk.Vertex(0, [1, 2, 3]),
            stk.Vertex(1, [-1, 2, -32]),
        ),
    ),
)
def vertices_(request) -> tuple[stk.Vertex, ...]:
    return request.param()
