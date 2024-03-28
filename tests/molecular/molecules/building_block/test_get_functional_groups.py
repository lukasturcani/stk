import itertools as it

from ..utilities import is_equivalent_functional_group

import stk
from collections.abc import Sequence
from .case_data import CaseData


def test_get_functional_groups(case_data: CaseData) -> None:
    """
    Test :meth:`.BuildingBlock.get_functional_groups`.

    Parameters:
        case_data:
            A test case. Holds the building block to test and the correct
            functional groups.

    """

    _test_get_functional_groups(
        building_block=case_data.building_block,
        functional_groups=case_data.functional_groups,
    )


def _test_get_functional_groups(
    building_block: stk.BuildingBlock,
    functional_groups: Sequence[stk.FunctionalGroup],
) -> None:
    """
    Test :meth:`.BuildingBlock.get_functional_groups`.

    Parameters:
        building_block:
            The building block to test.

        functional_groups:
            The correct functional groups.

    """

    for fg1, fg2 in it.zip_longest(
        building_block.get_functional_groups(),
        functional_groups,
    ):
        is_equivalent_functional_group(fg1, fg2)
