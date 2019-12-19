import pytest
import stk


@pytest.fixture(
    params=[
        stk.H(0),
    ],
)
def atom(request):
    return request.param.clone()
