import pytest
from pytest_lazy_fixtures import lf

from ...case_data import CaseData
from .helix import *

# Fixtures need to be visible for lazy_fixture() calls.
from .linear import *


@pytest.fixture(
    params=(
        lf("polymer_linear"),
        lf("polymer_helix"),
    ),
)
def polymer(request: pytest.FixtureRequest) -> CaseData:
    return request.param
