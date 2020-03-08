import pytest


@pytest.fixture(
    params=(
    ),
)
def cage_six_plus_nine(request):
    return request.param
