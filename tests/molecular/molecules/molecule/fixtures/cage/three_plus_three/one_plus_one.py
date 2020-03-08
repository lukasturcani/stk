import pytest


@pytest.fixture(
    params=(
    ),
)
def cage_one_plus_one(request):
    return request.param
