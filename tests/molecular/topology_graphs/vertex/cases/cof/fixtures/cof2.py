import pytest
import stk

from ....case_data import CaseData

vertices = stk.cof.honeycomb
vertices2 = stk.cof.vertices


@pytest.fixture
def cof2(init_at_center, id, vertices_, aligner_edge, cell):
    return CaseData(
        vertex=init_at_center(id, vertices_, aligner_edge, cell),
        id=id,
        position=(
            sum(v.get_position() for v in vertices_)/len(vertices_)
        ),
        cell=cell,
    )


@pytest.fixture(
    params=(
        vertices.LinearCofVertex.init_at_center,
        vertices.NonLinearCofVertex.init_at_center,
        vertices2.UnaligningVertex.init_at_center,
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
