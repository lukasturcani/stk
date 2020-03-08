import pytest


@pytest.fixture(
    params=(
    ),
)
def cage_two_plus_three(request):
    return request.param
