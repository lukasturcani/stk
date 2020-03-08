import pytest


@pytest.fixture(
    params=(
    ),
)
def cage_ten_plus_twenty(request):
    return request.param
