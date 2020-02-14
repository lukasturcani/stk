from ..utilities import has_same_structure, is_equivalent_molecule


def test_clone(molecule):
    clone = molecule.clone()
    has_same_structure(molecule, clone)
    assert molecule is not clone
    is_equivalent_molecule(molecule, clone)
