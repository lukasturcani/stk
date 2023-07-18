import pytest
from pytest_lazyfixture import lazy_fixture

from .case_data import CaseData

# Fixtures must be visible for lazy_fixture() calls.
from .fixtures import *  # noqa


@pytest.fixture(
    params=(
        lazy_fixture("add"),
        lazy_fixture("divide_by_mean"),
        lazy_fixture("multiply"),
        lazy_fixture("null"),
        lazy_fixture("power"),
        lazy_fixture("replace_fitness"),
        lazy_fixture("sequence"),
        lazy_fixture("shift_up"),
        lazy_fixture("sum"),
    ),
)
def case_data(request: pytest.FixtureRequest) -> CaseData:
    return request.param
