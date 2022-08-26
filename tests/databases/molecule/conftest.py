import pytest
from pytest_lazyfixture import lazy_fixture

# Fixtures need to be visible for lazy_fixture() calls.
from .fixtures import *  # noqa


@pytest.fixture(
    params=(lazy_fixture("molecule_mongo_db"),),
)
def case_data(request):
    return request.param
