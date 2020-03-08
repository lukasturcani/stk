import pytest


@pytest.fixture(
    params=(
    ),
)
def six_plus_twelve(request):
    return request.param
