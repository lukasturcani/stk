import pytest


@pytest.fixture(
    params=(
    ),
)
def cof_linkerless_honeycomb(request):
    return request.param
