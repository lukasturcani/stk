import pytest
from pytest_lazyfixture import lazy_fixture

# Fixtures must be visible for lazy_fixture() calls.
from .fixtures import *  # noqa


@pytest.fixture(
    params=(
        lazy_fixture('divide_by_mean'),
        lazy_fixture('multiply'),
        lazy_fixture('null'),
        lazy_fixture('power'),
    ),
)
def case_data(request):
    return request.param
