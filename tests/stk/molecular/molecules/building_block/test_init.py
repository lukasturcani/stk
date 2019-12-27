import rdkit.Chem.AllChem as rdkit
import stk

from .utilities import is_equivalent_building_block
from ..utilities import has_same_structure


def test_init(building_block, get_functional_groups):
    building_block = get_canonical_building_block(
        building_block=building_block,
        get_functional_groups=get_functional_groups,
    )
    rdkit_molecule = building_block.to_rdkit_mol()
    rdkit.SanitizeMol(rdkit_molecule)
    new = stk.BuildingBlock(
        smiles=rdkit.MolToSmiles(rdkit_molecule),
        functional_groups=building_block.get_functional_groups(),
    )
    is_equivalent_building_block(building_block, new)
    has_same_structure(building_block, new)


def get_canonical_building_block(
    building_block,
    get_functional_groups,
):
    rdkit_molecule = building_block.to_rdkit_mol()
    rdkit.SanitizeMol(rdkit_molecule)
    return stk.BuildingBlock(
        smiles=rdkit.MolToSmiles(rdkit_molecule, kekuleSmiles=True),
        functional_groups=get_functional_groups(building_block),
    )
