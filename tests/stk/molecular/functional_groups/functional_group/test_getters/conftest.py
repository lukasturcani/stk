import pytest

from .fixtures import *


@pytest.fixture(
    pytest.lazy_fixture('primary_amine'),
    pytest.lazy_fixture('secondary_amine'),
    pytest.lazy_fixture('aldehyde'),
    pytest.lazy_fixture('carboxylic_acid'),
)
def test_case(request):
    return request.param
