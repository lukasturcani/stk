import pytest


@pytest.fixture(
    params=(
    ),
)
def cage_two_plus_four(request):
    return request.param
