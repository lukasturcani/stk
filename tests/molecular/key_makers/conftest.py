import pytest
import stk

from .case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            key_maker=stk.Inchi(),
            molecule=stk.BuildingBlock('NCCN'),
            key_name='InChI',
            key='InChI=1S/C2H8N2/c3-1-2-4/h1-4H2',
        ),
    ),
)
def case_data(request):
    return request.param
