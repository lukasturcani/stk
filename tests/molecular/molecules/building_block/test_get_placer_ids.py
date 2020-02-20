import pytest
import itertools as it


def test_get_placer_ids(
    building_block,
    get_functional_groups,
):
    building_block = building_block.with_functional_groups(
        functional_groups=get_functional_groups(building_block),
    )
    placer_ids = it.zip_longest(
        building_block.get_placer_ids(),
        get_placer_ids(
            functional_groups=get_functional_groups(building_block),
        ),
    )
    for placer1, placer2 in placer_ids:
        assert placer1 == placer2


def get_placer_ids(functional_groups):
    for fg in functional_groups:
        yield from fg.get_placer_ids()
