import pytest


@pytest.fixture(
    params=(
    ),
)
def cage_five_plus_ten(request):
    return request.param
