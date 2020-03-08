import pytest


@pytest.fixutre(
    params=(
    ),
)
def cage_six_plus_eight(request):
    return request.param
