import pytest
from pytest_lazy_fixtures import lf

# Fixtures need to be visible for lazy_fixture() calls.
from .fixtures import *


@pytest.fixture(
    params=(
        lf("host"),
        lf("guest"),
    ),
)
def case_data(request):
    """
    A :class:`.CaseData` instance.

    """

    return request.param
