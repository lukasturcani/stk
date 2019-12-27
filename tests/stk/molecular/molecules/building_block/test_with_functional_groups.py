import itertools as it
import stk
import pytest

from ..utilities import has_same_structure
from .utilities import (
    is_equivalent_building_block,
    is_equivalent_functional_group,
)


class _TestCase:
    def __init__(self, building_block, factory, functional_groups):
        self.building_block = building_block
        self.factory = factory
        self.functional_groups = functional_groups


@pytest.mark.parametrize(
    argnames='test_case',
    argvalues=(
        _TestCase(
            building_block=stk.BuildingBlock(
                smiles='NCC(Br)CN',
                functional_groups=[stk.BromoFactory()],
            ),
            factory=stk.AmineFactory(),
            functional_groups=(
                stk.Amine(
                    atoms=(stk.N(0), stk.H(6), stk.H(7)),
                    bonders=(stk.N(0), ),
                    deleters=(stk.H(6), stk.H(7)),
                ),
                stk.Amine(
                    atoms=(stk.N(5), stk.H(13), stk.H(14)),
                    bonders=(stk.N(5), ),
                    deleters=(stk.H(13), stk.H(14)),
                ),
            ),
        ),
    ),
)
def test_with_functional_groups(test_case):
    _test_with_functional_groups(
        building_block=test_case.building_block,
        factory=test_case.factory,
        functional_groups=test_case.functional_groups,
    )


def _test_with_functional_groups(
    building_block,
    factory,
    functional_groups,
):
    clone = building_block.clone()
    _test_with_functional_groups_0(
        bulding_block=building_block,
        factory=factory,
        functional_groups=functional_groups,
    )
    # Test immutability.
    is_equivalent_building_block(building_block, clone)
    has_same_structure(building_block, clone)


def _test_with_functional_groups_0(
    building_block,
    factory,
    functional_groups
):
    new = building_block.with_functional_groups(factory)
    functional_groups = it.zip_longest(
        new.get_functional_groups(),
        it.chain(
            building_block.get_functional_groups(),
            functional_groups,
        ),
    )
    for fg1, fg2 in functional_groups:
        is_equivalent_functional_group(fg1, fg2)
