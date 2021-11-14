import stk

from .constructed_molecule import is_equivalent_constructed_molecule
from .molecule import is_equivalent_molecule


def is_equivalent(molecule1, molecule2):
    is_constructed_molecule1 = isinstance(
        molecule1,
        stk.ConstructedMolecule,
    )
    is_constructed_molecule2 = isinstance(
        molecule2,
        stk.ConstructedMolecule,
    )
    if is_constructed_molecule1 and is_constructed_molecule2:
        is_equivalent_constructed_molecule(molecule1, molecule2)
        return

    is_equivalent_molecule(molecule1, molecule2)
