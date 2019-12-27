import itertools as it
import stk
import pytest

from ..utilities import has_same_structure
from .utilities import (
    is_equivalent_building_block,
    is_equivalent_functional_group,
)


@pytest.fixture(
    params=(
        stk.Amine(
            atoms=(stk.N(0), stk.H(1), stk.H(7)),
            bonders=(stk.N(0), ),
            deleters=(stk.H(7), ),
        ),
        stk.Aldehyde(
            atoms=(stk.C(12), stk.O(2), stk.H(32)),
            bonders=(stk.C(12), ),
            deleters=(stk.H(32), stk.O(2)),
        ),
    ),
)
def functional_group(request):
    """
    A :class:`.FunctionalGroup` instance.

    """

    return request.param.clone()


def test_with_functional_group(building_block, functional_group):
    clone = building_block.clone()
    _test_with_functional_group(building_block, functional_group)
    # Test immutability.
    is_equivalent_building_block(building_block, clone)
    has_same_structure(building_block, clone)


def _test_with_functional_group(building_block, functional_group):
    new = building_block.with_functional_group(functional_group)
    original = it.chain(
        building_block.get_functional_groups(),
        (functional_group, ),
    )
    fgs = it.zip_longest(original, new.get_functional_groups)
    for original_fg, new_fg in fgs:
        is_equivalent_functional_group(original_fg, new_fg)
