from pytest_lazyfixture import lazy_fixture
import pytest

# Fixtures must be visible for lazy_fixture() calls.
from .nrotaxane import *  # noqa


@pytest.fixture(
    params=(
        lazy_fixture('nrotaxane'),
    ),
)
def rotaxane(request):
    return request.param
