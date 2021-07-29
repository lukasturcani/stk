import pytest
from pytest_lazyfixture import lazy_fixture

# Fixtures need to be visible for lazy_fixture() calls.
from .fixtures import *  # noqa
from .case_data import CaseData


@pytest.fixture(
    params=(
        lazy_fixture('molecule_mongo_dbs'),
    ),
)
def case_data(request) -> CaseData:
    return request.param