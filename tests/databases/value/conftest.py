import pytest
from pytest_lazyfixture import lazy_fixture

# Fixtures must be visible for lazy_fixture() calls.
from .fixtures import *  # noqa


@pytest.fixture(
    params=(
        lazy_fixture('constructed_molecule_value_mongo_db'),
        lazy_fixture('molecule_value_mongo_db'),
    ),
)
def case_data(request):
    return request.param
