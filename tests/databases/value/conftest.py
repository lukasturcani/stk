import pytest
from pytest_lazyfixture import lazy_fixture

# Fixtures must be visible for lazy_fixture() calls.
from .fixtures import *  # noqa


@pytest.fixture(
    params=(
        lazy_fixture('mongo_db_constructed_molecule_value_cache'),
        lazy_fixture('mongo_db_molecule_value_cache'),
        lazy_fixture('ram_constructed_molecule_value_cache'),
        lazy_fixture('ram_molecule_value_cache'),
    ),
)
def case_data(request):
    return request.param
