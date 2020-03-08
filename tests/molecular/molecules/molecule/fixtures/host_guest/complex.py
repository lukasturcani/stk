import pytest


@pytest.fixture(
    params=(
    ),
)
def host_guest_complex(request):
    return request.param
