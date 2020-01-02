from ..utilities import is_clone_molecule, has_same_structure


def test_clone(molecule):
    clone = molecule.clone()
    has_same_structure(molecule, clone)
    is_clone_molecule(molecule, clone)
