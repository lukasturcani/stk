import pytest
from pytest_lazyfixture import lazy_fixture

# Fixtures need to be visible for lazy_fixture() calls.
from .linear import *  # noqa


@pytest.fixture(
    params=(
        lazy_fixture('polymer_linear'),
        lazy_fixture('polymer_linear_opt'),
    ),
)
def polymer(request):
    return request.param
