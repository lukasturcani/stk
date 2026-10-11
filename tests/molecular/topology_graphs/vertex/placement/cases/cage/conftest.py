import pytest
from pytest_lazy_fixtures import lf

# Fixtures need to be visible for lazy_fixture() calls.
from .fixtures import *


@pytest.fixture(
    params=(
        lf("linear"),
        lf("nonlinear"),
        lf("unaligning"),
        lf("angled"),
    ),
)
def case_data(request):
    return request.param
