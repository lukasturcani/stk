from .case_data import CaseData


def test_get_num_functional_groups(case_data: CaseData) -> None:
    """
    Test :meth:`.BuildingBlock.get_num_functional_groups`.

    Parameters:
        case_data:
            A test case. Holds the building block to test and the correct
            functional groups.

    """

    assert case_data.building_block.get_num_functional_groups() == len(
        case_data.functional_groups
    )
