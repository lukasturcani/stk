import itertools as it

from .utilities import is_equivalent_functional_group


def test_get_functional_groups(building_block, get_functional_groups):
    building_block = building_block.with_functional_groups(
        functional_groups=get_functional_groups(building_block),
    )
    functional_groups = it.zip_longest(
        get_functional_groups(building_block),
        building_block.get_functional_groups(),
    )
    for fg1, fg2 in functional_groups:
        is_equivalent_functional_group(fg1, fg2)
