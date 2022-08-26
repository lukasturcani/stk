import pytest
from pytest_lazyfixture import lazy_fixture

# Fixtures need to be visible for lazy_fixture() calls.
from .fixtures import *  # noqa


@pytest.fixture(
    params=(
        lazy_fixture("linear"),
        lazy_fixture("nonlinear"),
        lazy_fixture("unaligning"),
    ),
)
def case_data(request):
    return request.param
