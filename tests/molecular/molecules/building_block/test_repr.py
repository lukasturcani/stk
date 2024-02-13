def test_repr(case_data):
    """
    Test :meth:`.BuildingBlock.__repr__`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the building block to test and the correct
        functional groups.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_repr(
        building_block=case_data.building_block,
        known_repr=case_data.known_repr,
    )


def _test_repr(building_block, known_repr):
    """
    Test :meth:`.BuildingBlock.get_functional_groups`.

    Parameters
    ----------
    building_block : :class:`.BuildingBlock`
        The building block to test.

    functional_groups : :class:`tuple` of :class:`.FunctionalGroup`
        The correct functional groups.

    Returns
    -------
    None : :class:`NoneType`

    """
    print(building_block.__repr__(), known_repr)
    assert building_block.__repr__() == known_repr
