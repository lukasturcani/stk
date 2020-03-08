import pytest


@pytest.fixture(
    params=(
    ),
)
def cof_hexagonal(request):
    return request.param
