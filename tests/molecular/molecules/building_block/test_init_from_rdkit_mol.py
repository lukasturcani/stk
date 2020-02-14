import stk

from .utilities import is_equivalent_building_block
from ..utilities import has_same_structure


def test_init_from_rdkit_mol(building_block):
    rdkit_molecule = building_block.to_rdkit_mol()
    new = stk.BuildingBlock.init_from_rdkit_mol(
        molecule=rdkit_molecule,
        functional_groups=building_block.get_functional_groups(),
    )
    is_equivalent_building_block(building_block, new)
    has_same_structure(building_block, new)
