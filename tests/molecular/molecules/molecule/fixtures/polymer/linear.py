import pytest


@pytest.fixture(
    params=(
    ),
)
def polymer_linear(request):
    return request.param
