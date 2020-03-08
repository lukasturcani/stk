import pytest
import pytest_lazyfixture as lazy_fixture

# Fixtures need to be visible for lazy_fixture() calls.
from .linear import *  # noqa


@pytest.fixture(
    params=(
        lazy_fixture('polymer_linear'),
    ),
)
def polymer(request):
    return request.param
