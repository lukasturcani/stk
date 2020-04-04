import pytest
from pytest_lazyfixture import lazy_fixture

# Fixtures must be visible for lazy_fixture() calls.
from .fixtures import *  # noqa


@pytest.fixture(
    params=(
        lazy_fixture('divide_by_mean'),
        lazy_fixture('multiply'),
        lazy_fixture('null'),
        lazy_fixture('replace_fitness'),
        lazy_fixture('shift_up'),
        lazy_fixture('sum'),
    ),
)
def case_data(request):
    return request.param
