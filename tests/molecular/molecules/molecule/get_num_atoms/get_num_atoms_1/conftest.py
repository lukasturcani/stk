import pytest
import stk

from .case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            molecule=stk.BuildingBlock('NCCN'),
            num_atoms=12,
        ),
    ),
)
def case_data(request):
    return request.param
