import pytest
from pytest_lazy_fixtures import lf

# Fixtures must be visible for lazy_fixture() calls.
from .fixtures import *


@pytest.fixture(
    params=(
        lf("genetic_recombination"),
        lf("random_crosser"),
    ),
)
def case_data(request):
    return request.param
