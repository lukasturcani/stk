import pytest
import stk

from ....case_data import CaseData


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


@pytest.fixture(
    params=(
        stk.cof.LinearVertex,
        stk.cof.NonLinearVertex,
        stk.cof.UnaligningVertex,
    ),
)
def cls(request):
    return request.param
