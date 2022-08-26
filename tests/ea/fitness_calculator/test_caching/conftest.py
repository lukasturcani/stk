import pytest
from pytest_lazyfixture import lazy_fixture

# Fixtures must be visible for lazy_fixture() calls.
from .fixtures import *  # noqa


@pytest.fixture(
    params=(
        lazy_fixture("fitness_function"),
        lazy_fixture("property_vector"),
    ),
)
def case_data(request):
    return request.param
