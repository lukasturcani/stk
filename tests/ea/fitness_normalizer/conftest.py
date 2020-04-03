import pytest
from pytest_lazyfixture import lazy_fixture


@pytest.fixture(
    params=(
        lazy_fixture('divide_by_mean'),
    ),
)
def case_data(request):
    return request.param
