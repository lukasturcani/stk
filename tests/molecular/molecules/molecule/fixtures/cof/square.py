import pytest


@pytest.fixture(
    params=(
    ),
)
def cof_square(request):
    return request.param
