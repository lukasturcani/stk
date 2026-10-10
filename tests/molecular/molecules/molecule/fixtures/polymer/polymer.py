import pytest
from pytest_lazyfixture import lazy_fixture

from ...case_data import CaseData
from .helix import *

# Fixtures need to be visible for lazy_fixture() calls.
from .linear import *


@pytest.fixture(
    params=(
        lazy_fixture("polymer_linear"),
        lazy_fixture("polymer_helix"),
    ),
)
def polymer(request: pytest.FixtureRequest) -> CaseData:
    return request.param
