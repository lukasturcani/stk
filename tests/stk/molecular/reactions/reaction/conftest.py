import pytest

from .fixtures import *


@pytest.fixture(
    params=(
        pytest.lazy_fixture('one_one_reaction'),
        pytest.lazy_fixture('one_two_reaction'),
        pytest.lazy_fixture('two_two_reaction'),
        pytest.lazy_fixture('ring_amine_reaction'),
    ),
)
def test_case(request):
    return request.param
