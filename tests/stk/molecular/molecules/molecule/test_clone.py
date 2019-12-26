from .utilities import is_equivalent_molecule, has_same_structure


def test_clone(molecule):
    clone = molecule.clone()
    is_equivalent_molecule(molecule, clone)
    has_same_structure(molecule, clone)
