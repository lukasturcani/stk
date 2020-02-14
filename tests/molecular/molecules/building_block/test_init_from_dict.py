from .utilities import is_equivalent_building_block
from ..utilities import has_same_structure


def test_init_from_dict(building_block):
    new = building_block.__class__.init_from_dict(
        molecule_dict=building_block.to_dict(),
    )
    is_equivalent_building_block(building_block, new)
    has_same_structure(building_block, new)
