import pytest
import stk

from ....case_data import CaseData

vertices = stk.framework.honeycomb


@pytest.fixture
def framework1(cls, id, position, aligner_edge, cell):
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
        vertices._LinearFrameworkVertex,
        vertices._NonLinearFrameworkVertex,
    ),
)
def cls(request):
    return request.param
