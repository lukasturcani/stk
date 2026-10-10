import pytest
from pytest_lazy_fixtures import lf

# Fixtures must be visible for lazy_fixture() calls.
from .nrotaxane import *


@pytest.fixture(
    params=(lf("rotaxane_nrotaxane"),),
)
def rotaxane(request):
    return request.param
