import pytest
from pytest_lazy_fixtures import lf

from .case_data import CaseData

# Fixtures need to be visible for lazy_fixture() calls.
from .fixtures import *


@pytest.fixture(
    params=(lf("constructed_molecule_mongo_db"),),
)
def case_data(request) -> CaseData:
    return request.param
