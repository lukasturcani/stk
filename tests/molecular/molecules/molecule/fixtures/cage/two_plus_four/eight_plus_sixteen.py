import pytest


@pytest.fixture(
    params=(
    ),
)
def cage_eight_plus_sixteen(request):
    return request.param
