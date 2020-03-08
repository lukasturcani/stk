import pytest


@pytest.fixture(
    params=(
    ),
)
def cage_four_plus_six(request):
    return request.param
