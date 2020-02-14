import stk

from .utilities import is_equivalent_building_block
from ..utilities import has_same_structure


def test_load(building_block, tmpdir):
    path = str(tmpdir / 'building_block.json')
    building_block.dump(path)
    new = stk.BuildingBlock.load(path)
    is_equivalent_building_block(building_block, new)
    has_same_structure(building_block, new)
