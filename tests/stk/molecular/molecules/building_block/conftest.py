import stk
import pytest


@pytest.fixture(
    params=[
        stk.BuildingBlock('NCCN'),
        stk.BuildingBlock('NCCN', [stk.AmineFactory()]),
        stk.BuildingBlock(
            smiles='BrCC(Br)C(Br)C(Br)C(Br)C(Br)C(Br)CBr',
            functional_groups=[stk.BromoFactory()],
        ),
        stk.BuildingBlock('N[C+][C+2]N'),
    ],
)
def building_block(request):
    return request.param.clone()


@pytest.fixture(
    params=[
        lambda building_block: None,
        lambda building_block: range(
            building_block.get_num_functional_groups()
        ),
        lambda building_block:
            range(0, building_block.get_num_functional_groups(), 2),
        lambda building_block:
            range(min(1, building_block.get_num_functional_groups())),
        lambda building_block:
            list(range(min(
                1, building_block.get_num_functional_groups()
            ))),
        lambda building_block:
            tuple(range(min(
                1, building_block.get_num_functional_groups()
            ))),
        lambda building_block: (
            i for i in range(min(
                1, building_block.get_num_functional_groups()
            ))
        ),
    ],
)
def get_fg_ids(request):
    return request.param
