import itertools as it
from .utilities import is_atom_clone


def test_get_atoms(test_case):
    _test_get_atoms(
        functional_group=test_case.functional_group,
        atoms=test_case.atoms,
    )


def _test_get_atoms(functional_group, atoms):
    fg_atoms = it.zip_longest(functional_group.get_atoms(), atoms)
    for atom1, atom2 in fg_atoms:
        is_atom_clone(atom1, atom2)
