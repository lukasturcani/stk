from tests.utilities import is_equivalent_constructed_molecule


def is_clone_constructed_molecule(
    constructed_molecule1,
    constructed_molecule2,
):
    assert constructed_molecule1 is not constructed_molecule2
    is_equivalent_constructed_molecule(
        constructed_molecule1=constructed_molecule1,
        constructed_molecule2=constructed_molecule2,
    )
