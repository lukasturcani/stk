import pytest


@pytest.fixture(
    params=(
    ),
)
def cage_four_plus_eight(request):
    return request.param
