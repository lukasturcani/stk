import itertools as it
from .utilities import is_atom_clone


def test_get_bonders(test_case):
    _test_get_bonders(
        functional_group=test_case.functional_group,
        bonders=test_case.bonders,
    )


def _test_get_bonders(functional_group, bonders):
    fg_atoms = it.zip_longest(functional_group.get_bonders(), bonders)
    for atom1, atom2 in fg_atoms:
        is_atom_clone(atom1, atom2)
