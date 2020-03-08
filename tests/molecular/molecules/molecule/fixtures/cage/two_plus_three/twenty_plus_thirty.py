import pytest


@pytest.fixture(
    params=(
    ),
)
def cage_twenty_plus_thirty(request):
    return request.param
