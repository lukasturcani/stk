from ..utilities import (
    are_equivalent_functional_groups,
    has_same_structure,
    is_equivalent_building_block,
    is_equivalent_molecule,
)


def test_with_functional_groups(building_block, get_functional_groups):
    """
    Test :meth:`.BuildingBlock.with_functional_groups`.

    Parameters
    ----------
    building_block : :class:`.BuildingBlock`
        The building block to test.

    get_functional_groups : :class:`callable`
        Takes a single parameter, `building_block` and returns the
        `functional_groups` parameter to use for this test.

    Returns
    -------
    None : :class:`NoneType`

    """

    # Save clone to check immutability.
    clone = building_block.clone()
    _test_with_functional_groups(
        building_block=building_block,
        functional_groups=tuple(get_functional_groups(building_block)),
    )
    is_equivalent_building_block(building_block, clone)
    has_same_structure(building_block, clone)


def _test_with_functional_groups(building_block, functional_groups):
    """
    Test :meth:`.BuildingBlock.with_functional_groups`.

    Parameters
    ----------
    building_block : :class:`.BuildingBlock`
        The building block to test.

    functional_groups : :class:`tuple` of :class:`.FunctionalGroup`
        The functional groups the new building block should hold.

    Returns
    -------
    None : :class:`NoneType`

    """

    new = building_block.with_functional_groups(functional_groups)
    are_equivalent_functional_groups(
        new.get_functional_groups(),
        functional_groups,
    )
    is_equivalent_molecule(building_block, new)
    has_same_structure(building_block, new)
