import pytest
import stk


from ...case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            molecule=stk.BuildingBlock('NCCN'),
            smiles='NCCN',
        ),
    ),
)
def case_data(request):
    return request.param
