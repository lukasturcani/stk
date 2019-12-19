import pytest
import stk


@pytest.fixture(
    params=[
        stk.Bond(stk.H(0), stk.H(1), 1),
        stk.Bond(stk.He(2), stk.He(3), 2, (1, 0, -1)),
    ],
)
def bond(request):
    return request.param.clone()
