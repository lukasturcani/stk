def test_get_num_functional_groups(case_data):
    """
    Test :meth:`.BuildingBlock.get_num_functional_groups`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the building block to test and the
        correct number of functional groups.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert case_data.building_block.get_num_functional_groups() == len(
        case_data.functional_groups
    )
