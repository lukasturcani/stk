def test_get_building_block_id(case_data):
    """
    Test :meth:`.AtomInfo.get_building_block_id`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the atom info to test and the correct
        building block id.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert (
        case_data.atom_info.get_building_block_id()
        == case_data.building_block_id
    )
