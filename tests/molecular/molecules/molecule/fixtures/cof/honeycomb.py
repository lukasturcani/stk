import pytest


@pytest.fixture(
    params=(
    ),
)
def cof_honeycomb(request):
    return request.param
