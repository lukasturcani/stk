import pytest
from pytest_lazy_fixtures import lf

# Fixtures need to be visible for lazy_fixture() calls.
from .fixtures import *


@pytest.fixture(
    params=(
        lf("linear"),
        lf("cof"),
    ),
)
def case_data(request):
    return request.param


@pytest.fixture
def constructed_molecule(case_data):
    return case_data.constructed_molecule
