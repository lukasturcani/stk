import pytest
from pytest_lazy_fixtures import lf

# Fixtures need to be visible for lazy_fixture() calls.
from .complex import *


@pytest.fixture(
    params=(lf("host_guest_complex"),),
)
def host_guest(request):
    return request.param
