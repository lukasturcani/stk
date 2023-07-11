import pytest
from pytest_lazyfixture import lazy_fixture

from .case_data import CaseData

# Fixtures must be visible for lazy_fixture() calls.
from .fixtures import *  # noqa


@pytest.fixture(
    params=(
        lazy_fixture("above_average"),
        lazy_fixture("best"),
        lazy_fixture("filter_batches"),
        lazy_fixture("filter_molecules"),
        lazy_fixture("remove_batches"),
        lazy_fixture("remove_molecules"),
        lazy_fixture("worst"),
    ),
)
def case_data(request: pytest.FixtureRequest) -> CaseData:
    return request.param
