import pytest
from pytest_lazyfixture import lazy_fixture

# Fixtures need to be visible for lazy_fixture() calls.
from .fixtures import *  # noqa


@pytest.fixture(
    params=(
        lazy_fixture("host"),
        lazy_fixture("guest"),
    ),
)
def case_data(request):
    """
    A :class:`.CaseData` instance.

    """

    return request.param
