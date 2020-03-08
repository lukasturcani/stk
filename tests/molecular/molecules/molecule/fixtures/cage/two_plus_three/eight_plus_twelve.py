import pytest


@pytest.fixture(
    params=(
    ),
)
def cage_eight_plus_twelve(request):
    return request.param
