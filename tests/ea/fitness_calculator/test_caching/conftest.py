import pytest
from pytest_lazy_fixtures import lf

from .case_data import CaseData

# Fixtures must be visible for lazy_fixture() calls.
from .fixtures import *


@pytest.fixture(
    params=(
        lf("fitness_function"),
        lf("property_vector"),
    ),
)
def case_data(request: pytest.FixtureRequest) -> CaseData:
    return request.param
