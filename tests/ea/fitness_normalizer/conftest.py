import pytest
from pytest_lazy_fixtures import lf

from .case_data import CaseData

# Fixtures must be visible for lazy_fixture() calls.
from .fixtures import *


@pytest.fixture(
    params=(
        lf("add"),
        lf("divide_by_mean"),
        lf("multiply"),
        lf("null"),
        lf("power"),
        lf("replace_fitness"),
        lf("sequence"),
        lf("shift_up"),
        lf("sum"),
    ),
)
def case_data(request: pytest.FixtureRequest) -> CaseData:
    return request.param
