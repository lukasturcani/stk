import stk
import itertools as it
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
    params=(
        lambda molecule:
            stk.BromoFactory().get_functional_groups(molecule),
        lambda molecule:
            stk.AmineFactory().get_functional_groups(molecule),
        lambda molecule:
            it.chain(
                stk.AmineFactory().get_functional_groups(molecule),
                stk.BromoFactory().get_functional_groups(molecule),
            )
    )
)
def get_functional_groups(request):
    return request.param
