from .case_data import CaseData

import stk


def test_repr(case_data: CaseData) -> None:
    """
    Test :meth:`.BuildingBlock.__repr__`.

    Parameters:
        case_data:
            A test case. Holds the building block to test.

    """

    _test_repr(
        building_block=case_data.building_block,
        known_repr=case_data.known_repr,
    )


def _test_repr(building_block: stk.BuildingBlock, known_repr: str) -> None:
    """
    Test :meth:`.BuildingBlock.get_functional_groups`.

    Parameters:
        building_block:
            The building block to test.

        known_repr:
            The correct representation.

    """

    print(building_block.__repr__(), known_repr)
    assert building_block.__repr__() == known_repr
