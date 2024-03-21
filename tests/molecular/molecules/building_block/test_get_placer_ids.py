import itertools as it
from .case_data import CaseData

import stk


def test_get_placer_ids(case_data: CaseData) -> None:
    """
    Test :meth:`.BuildingBlock.get_placer_ids`.

    Parameters:
        case_data:
            A test case. Holds the building block to test and the
            correct placer ids.

    """

    _test_get_placer_ids(
        building_block=case_data.building_block,
        placer_ids=case_data.placer_ids,
    )


def _test_get_placer_ids(
    building_block: stk.BuildingBlock,
    placer_ids: tuple[int, ...],
) -> None:
    """
    Test :meth:`.BuildingBlock.get_placer_ids`.

    Parameters:
        building_block:
            The building block to test.

        placer_ids:
            The correct placer ids.

    """

    for placer1, placer2 in it.zip_longest(
        building_block.get_placer_ids(),
        placer_ids,
    ):
        assert placer1 == placer2
