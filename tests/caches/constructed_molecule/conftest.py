import pytest
from pytest_lazyfixture import lazy_fixture

# Fixtures need to be visible for lazy_fixture() calls.
from .fixtures import *  # noqa


@pytest.fixture(
    params=(
        lazy_fixture('mongo_db_constructed_molecule_cache'),
        lazy_fixture('ram_constructed_molecule_cache'),
    ),
)
def case_data(request):
    return request.param
