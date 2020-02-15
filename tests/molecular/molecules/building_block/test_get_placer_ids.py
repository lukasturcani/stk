import pytest
import itertools as it


@pytest.fixture(
    params=[
        lambda building_block: None,
        lambda building_block: (
            0
            if building_block.get_num_functional_groups() > 0
            else None
        ),
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


def test_get_placer_ids(
    building_block,
    get_functional_groups,
    get_fg_ids,
):
    building_block = building_block.with_functional_groups(
        functional_groups=get_functional_groups(building_block),
    )
    placer_ids = it.zip_longest(
        building_block.get_placer_ids(get_fg_ids(building_block)),
        get_placer_ids(
            functional_groups=get_functional_groups(building_block),
            fg_ids=normalize_ids(
                building_block=building_block,
                ids=get_fg_ids(building_block),
            ),
        ),
    )
    for placer1, placer2 in placer_ids:
        assert placer1 == placer2


def get_placer_ids(functional_groups, fg_ids):
    fg_ids = set(fg_ids)
    for fg_id, fg in enumerate(functional_groups):
        if fg_id in fg_ids:
            yield from fg.get_placer_ids()


def normalize_ids(building_block, ids):
    if ids is None:
        return range(building_block.get_num_functional_groups())
    elif isinstance(ids, int):
        return (ids, )
    else:
        return ids
