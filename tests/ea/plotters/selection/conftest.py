import pytest

from .case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
        ),
    ),
)
def case_data(request):
    return request.param
