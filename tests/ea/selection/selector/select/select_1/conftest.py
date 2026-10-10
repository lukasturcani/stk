import pytest
from pytest_lazy_fixtures import lf

from .case_data import CaseData

# Fixtures must be visible for lazy_fixture() calls.
from .fixtures import *


@pytest.fixture(
    params=(
        lf("above_average"),
        lf("best"),
        lf("filter_batches"),
        lf("filter_molecules"),
        lf("remove_batches"),
        lf("remove_molecules"),
        lf("worst"),
    ),
)
def case_data(request: pytest.FixtureRequest) -> CaseData:
    return request.param
