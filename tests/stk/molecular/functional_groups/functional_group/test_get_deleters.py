import itertools as it

from .utilities import is_atom_clone


def test_get_deleters(test_case):
    _test_get_deleters(
        functional_group=test_case.functional_group,
        deleters=test_case.deleters,
    )


def _test_get_deleters(functional_group, deleters):
    fg_atoms = it.zip_longest(
        functional_group.get_deleters(),
        deleters,
    )
    for atom1, atom2 in fg_atoms:
        is_atom_clone(atom1, atom2)
