import pytest


@pytest.fixture(
    params=(
    ),
)
def cage_two_plus_two(request):
    return request.param
