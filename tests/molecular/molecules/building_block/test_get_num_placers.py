import stk

from .case_data import CaseData


def test_get_num_placers(case_data: CaseData) -> None:
    """
    Test :meth:`.BuildingBlock.get_num_placers`.

    Parameters:

        case_data:
            A test case. Holds the building block to test and the
            correct number of *placer* atoms.

    """

    _test_get_num_placers(
        building_block=case_data.building_block,
        expected_num_placers=len(case_data.placer_ids),
    )


def _test_get_num_placers(
    building_block: stk.BuildingBlock,
    expected_num_placers: int,
) -> None:
    """
    Test :meth:`.BuildingBlock.get_num_placers`.

    Parameters:

        building_block:
            The building block to test.

        expected_num_placers:
            The correct number of *placer* atoms

    """

    assert building_block.get_num_placers() == expected_num_placers
