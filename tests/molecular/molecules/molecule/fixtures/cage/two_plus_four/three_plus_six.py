import pytest


@pytest.fixture(
    params=(
    ),
)
def cage_three_plus_six(request):
    return request.param
