from .utilities import is_equivalent_building_block
from ..utilities import has_same_structure


def test_clone(building_block):
    clone = building_block.clone()
    is_equivalent_building_block(building_block, clone)
    has_same_structure(building_block, clone)
