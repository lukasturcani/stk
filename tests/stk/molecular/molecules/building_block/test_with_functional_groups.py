import stk
import pytest

from ..utilities import has_same_structure, is_equivalent_molecule
from .utilities import (
    is_equivalent_building_block,
    are_equivalent_functional_groups,
)


class _TestCase:
    def __init__(self, building_block, functional_groups):
        self.building_block = building_block
        self.functional_groups = functional_groups


@pytest.mark.parametrize(
    argnames='test_case',
    argvalues=(
        _TestCase(
            building_block=stk.BuildingBlock(
                smiles='NCC(Br)CN',
                functional_groups=[stk.BromoFactory()],
            ),
            functional_groups=(
                stk.PrimaryAmino(
                    nitrogen=stk.N(0),
                    hydrogen1=stk.H(6),
                    hydrogen2=stk.H(7),
                    atom=stk.C(1),
                    bonders=(stk.N(0), ),
                    deleters=(stk.H(6), stk.H(7)),
                ),
                stk.PrimaryAmino(
                    nitrogen=stk.N(5),
                    hydrogen1=stk.H(13),
                    hydrogen2=stk.H(14),
                    atom=stk.C(4),
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
        functional_groups=test_case.functional_groups,
    )


def _test_with_functional_groups(building_block, functional_groups):
    clone = building_block.clone()
    _test_with_functional_groups_0(
        building_block=building_block,
        functional_groups=functional_groups,
    )
    # Test immutability.
    is_equivalent_building_block(building_block, clone)
    has_same_structure(building_block, clone)


def _test_with_functional_groups_0(building_block, functional_groups):
    new = building_block.with_functional_groups(functional_groups)
    are_equivalent_functional_groups(
        new.get_functional_groups(),
        functional_groups,
    )
    is_equivalent_molecule(building_block, new)
    has_same_structure(building_block, new)
