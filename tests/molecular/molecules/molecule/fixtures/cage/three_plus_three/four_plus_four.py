import pytest


@pytest.fixture(
    params=(
    ),
)
def cage_four_plus_four(request):
    return request.param
