from .utilities import is_equivalent_building_block
from ..utilities import has_same_structure


def test_clone(building_block):
    """
    Test :meth:`.BuildingBlock.clone`.

    Parameters
    ----------
    building_block : :class:`.BuildingBlock`
        The building block to test.

    Returns
    -------
    None : :class:`NoneType`

    """

    clone = building_block.clone()
    assert building_block is not clone
    is_equivalent_building_block(building_block, clone)
    has_same_structure(building_block, clone)
