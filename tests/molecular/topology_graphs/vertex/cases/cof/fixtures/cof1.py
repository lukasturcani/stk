import pytest
import stk

from ....case_data import CaseData

vertices = stk.cof.honeycomb
vertices2 = stk.cof.vertices


@pytest.fixture
def cof1(cls, id, position, aligner_edge, cell):
    return CaseData(
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


@pytest.fixture
def unaligned_cof1(cls, id, position, aligner_edge, cell):
    return CaseData(
        vertex=vertices2._UnaligningVertex(
            vertex=cls(
                id=id,
                position=position,
                aligner_edge=aligner_edge,
                cell=cell,
            )
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
