import pytest


@pytest.fixture(
    params=(
    ),
)
def test_case(request):
    return request.param
