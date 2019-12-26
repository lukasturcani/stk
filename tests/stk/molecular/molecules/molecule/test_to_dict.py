import stk

from .utilities import has_same_structure, is_equivalent_molecule


def test_to_dict(molecule):
    new = stk.BuildingBlock.init_from_dict(molecule.to_dict())
    is_equivalent_molecule(molecule, new)
    has_same_structure(molecule, new)
