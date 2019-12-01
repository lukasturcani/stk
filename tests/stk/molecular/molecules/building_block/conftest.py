import stk
import pytest


@pytest.fixture(
    params=[
        stk.BuildingBlock('NCCN'),
        stk.BuildingBlock('NCCN', ['amine']),
        stk.BuildingBlock('N[C+][C+2]N'),
    ]
)
def building_block(request):
    return request.param.clone()
