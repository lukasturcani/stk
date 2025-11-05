import pytest
from pytest_lazyfixture import lazy_fixture

from ...case_data import CaseData

# Fixtures need to be visible for lazy_fixture() calls.
from .internal_reaction import *  # noqa
from .ncore import *  # noqa


@pytest.fixture(
    params=(
        lazy_fixture("small_ncore"),
        lazy_fixture("small_internal"),
    ),
)
def small(request: pytest.FixtureRequest) -> CaseData:
    return request.param
