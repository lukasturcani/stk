import pytest


@pytest.fixture(
    params=(
    ),
)
def cage_four_plus_six_2(request):
    return request.param
