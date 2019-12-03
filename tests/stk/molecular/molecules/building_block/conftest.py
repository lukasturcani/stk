import stk
import pytest


@pytest.fixture(
    params=[
        stk.BuildingBlock('NCCN'),
        stk.BuildingBlock('NCCN', ['amine']),
        stk.BuildingBlock('N[C+][C+2]N'),
    ],
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
    ],
)
def molecule(request):
    return request.param.clone()


@pytest.fixture(
    params=[
        lambda building_block: None,
        lambda building_block: range(len(building_block.func_groups)),
        lambda building_block:
            range(0, len(building_block.func_groups), 2),
        lambda building_block:
            range(0, min(1, len(building_block.func_groups))),
        lambda building_block:
            list(range(0, min(1, len(building_block.func_groups)))),
        lambda building_block:
            tuple(range(0, min(1, len(building_block.func_groups)))),
        lambda building_block: (
            i
            for i in range(0, min(1, len(building_block.func_groups)))
        ),
    ],
)
def get_fg_ids(request):
    return request.param
