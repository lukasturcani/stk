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


@pytest.fixture(
    params=[
        stk.BuildingBlock('NCCN'),
        stk.BuildingBlock('N[C+][C+2]N'),
        stk.ConstructedMolecule(
            building_blocks=[stk.BuildingBlock('BrCCBr', ['bromine'])],
            topology_graph=stk.polymer.Linear('A', 3),
        ),
    ]
)
def molecule(request):
    return request.param.clone()
