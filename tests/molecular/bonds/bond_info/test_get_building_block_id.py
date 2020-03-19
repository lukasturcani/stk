def test_get_building_block_id(case_data):
    """
    Test :meth:`.BondInfo.get_building_block_id`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the bond info to test and the building
        block id it should be holding.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert (
        case_data.bond_info.get_building_block_id()
        == case_data.building_block_id
    )
