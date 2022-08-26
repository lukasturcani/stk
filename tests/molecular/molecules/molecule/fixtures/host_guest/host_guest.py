import pytest
from pytest_lazyfixture import lazy_fixture

# Fixtures need to be visible for lazy_fixture() calls.
from .complex import *  # noqa


@pytest.fixture(
    params=(lazy_fixture("host_guest_complex"),),
)
def host_guest(request):
    return request.param
