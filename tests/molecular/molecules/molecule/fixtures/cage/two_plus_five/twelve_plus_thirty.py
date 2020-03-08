import pytest


@pytest.fixture(
    params=(
    ),
)
def cage_twelve_plus_thirty(request):
    return request.param
