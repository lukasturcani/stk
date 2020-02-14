import stk
import itertools as it
import pytest

from .utilities import is_equivalent_functional_group
from ..utilities import is_equivalent_molecule, has_same_structure


@pytest.fixture(
    params=(
        stk.BuildingBlock('NCCN'),
        stk.BuildingBlock('NCCN', [stk.BromoFactory()]),
        stk.BuildingBlock('BrCC(Br)CC(Br)CCBr'),
        stk.BuildingBlock('NCCBr'),
    )
)
def molecule(request):
    return request.param.clone()


def test_init_from_molecule(molecule, get_functional_groups):
    building_block = stk.BuildingBlock.init_from_molecule(
        molecule=molecule,
        functional_groups=get_functional_groups(molecule),
    )
    has_same_structure(molecule, building_block)
    is_equivalent_molecule(building_block, molecule)
    functional_groups = it.zip_longest(
        building_block.get_functional_groups(),
        get_functional_groups(molecule),
    )
    for fg1, fg2 in functional_groups:
        is_equivalent_functional_group(fg1, fg2)
