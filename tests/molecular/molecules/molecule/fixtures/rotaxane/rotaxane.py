from pytest_lazyfixture import lazy_fixture
import pytest


@pytest.fixture(
    params=(
        lazy_fixture('nrotaxane'),
    ),
)
def rotaxane(request):
    return request.param
