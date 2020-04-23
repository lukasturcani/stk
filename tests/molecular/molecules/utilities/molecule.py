from tests.utilities import is_equivalent_molecule


def is_clone_molecule(molecule1, molecule2):
    assert molecule1 is not molecule2
    is_equivalent_molecule(molecule1, molecule2)
