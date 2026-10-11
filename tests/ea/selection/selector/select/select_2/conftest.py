import pytest
from pytest_lazy_fixtures import lf

from .case_data import CaseData

# Fixtures must be visible for lazy_fixture() calls.
from .fixtures import *


@pytest.fixture(
    params=(
        lf("roulette"),
        lf("stochastic_universal_sampling"),
        lf("tournament"),
    ),
)
def case_data(request: pytest.FixtureRequest) -> CaseData:
    return request.param
