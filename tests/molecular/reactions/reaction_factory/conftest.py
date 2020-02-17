import pytest
from pytest_lazyfixture import lazy_fixture

# Fixtures must be visible for lazy_fixture() calls.
from .fixtures import *  # noqa


@pytest.fixture(
    params=(
        lazy_fixture('one_one_reaction'),
        lazy_fixture('one_two_reaction'),
        lazy_fixture('two_two_reaction'),
    ),
)
def test_case(request):
    return request.param
