import pytest
import stk


@pytest.mark.parametrize(
    argnames=('building_block', 'num_functional_groups'),
    argvalues=(
        (stk.BuildingBlock('NCCN'), 0),
        (stk.BuildingBlock('NCCN', [stk.PrimaryAmineFactory()]), 2),
        (stk.BuildingBlock('NCCN', [stk.BromoFactory()]), 0),
    )
)
def test_get_num_functional_groups(
    building_block,
    num_functional_groups,
):
    assert (
        building_block.get_num_functional_groups()
        == num_functional_groups
    )
