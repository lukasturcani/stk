import itertools as it


def test_get_placer_ids(case_data):
    """
    Test :meth:`.BuildingBlock.get_placer_ids`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the building block to test and the
        correct placer ids.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_get_placer_ids(
        building_block=case_data.building_block,
        placer_ids=case_data.placer_ids,
    )


def _test_get_placer_ids(building_block, placer_ids):
    """
    Test :meth:`.BuildingBlock.get_placer_ids`.

    Parameters
    ----------
    building_block : :class:`.BuildingBlock`
        The building block to test.

    placer_ids : :class:`tuple` of :class:`int`
        The correct placer ids.

    Returns
    -------
    None : :class:`NoneType`

    """

    for placer1, placer2 in it.zip_longest(
        building_block.get_placer_ids(),
        placer_ids,
    ):
        assert placer1 == placer2
