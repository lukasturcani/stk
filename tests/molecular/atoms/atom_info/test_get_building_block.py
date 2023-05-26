def test_get_building_block(case_data):
    """
    Test :meth:`.AtomInfo.get_building_block`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the atom info to test and the building block
        it should be holding.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert case_data.atom_info.get_building_block() is case_data.building_block
