import pytest
from pytest_lazy_fixtures import lf

from ...case_data import CaseData

# Fixtures need to be visible for lazy_fixture() calls.
from .internal_reaction import *
from .ncore import *


@pytest.fixture(
    params=(
        lf("small_ncore"),
        lf("small_internal"),
    ),
)
def small(request: pytest.FixtureRequest) -> CaseData:
    return request.param
