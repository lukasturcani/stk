import pytest
from pytest_lazyfixture import lazy_fixture

# Fixtures must be visible for lazy_fixture() calls.
from .nrotaxane import *  # noqa


@pytest.fixture(
    params=(lazy_fixture("rotaxane_nrotaxane"),),
)
def rotaxane(request):
    return request.param
