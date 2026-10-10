import pytest
from pytest_lazy_fixtures import lf

# Fixtures must be visible for lazy_fixture() calls.
from .fixtures import *


@pytest.fixture(
    params=(
        lf("random_building_block"),
        lf("random_topology_graph"),
        lf("similar_building_block"),
        lf("random_mutator"),
    ),
)
def case_data(request):
    return request.param
