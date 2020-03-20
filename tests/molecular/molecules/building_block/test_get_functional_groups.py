import itertools as it

from ..utilities import is_equivalent_functional_group


def test_get_functional_groups(case_data):
    """
    Test :meth:`.BuildingBlock.get_functional_groups`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the building block to test and the correct
        functional groups.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_get_functional_groups(
        building_block=case_data.building_block,
        functional_groups=case_data.functional_groups,
    )


def _test_get_functional_groups(building_block, functional_groups):
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

    for fg1, fg2 in it.zip_longest(
        building_block.get_functional_groups(),
        functional_groups,
    ):
        is_equivalent_functional_group(fg1, fg2)
