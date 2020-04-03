import pytest
import stk

from ..case_data import CaseData


@pytest.fixture(
    params=(
    ),
)
def property_vector(request):
    return request.param
