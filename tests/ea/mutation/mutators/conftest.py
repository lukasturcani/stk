import pytest
from pytest_lazyfixture import lazy_fixture


@pytest.fixture(
    params=(
        lazy_fixture('random_building_block'),
    ),
)
def case_data(request):
    return request.param
