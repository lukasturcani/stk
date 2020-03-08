import pytest


@pytest.fixture(
    params=(
    ),
)
def cof_kagome(request):
    return request.param
