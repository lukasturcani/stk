import pytest
import stk


from ...case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            molecule=stk.BuildingBlock('NCCN'),
            smiles='NCCN',
        ),
        CaseData(
            molecule=stk.BuildingBlock('[H]NCCN'),
            smiles='NCCN',
        ),
        CaseData(
            molecule=stk.BuildingBlock('C(N)CN'),
            smiles='NCCN',
        ),
    ),
)
def case_data(request):
    return request.param
