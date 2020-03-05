import itertools as it


def test_get_atoms(case_data):
    _test_get_atoms(
        functional_group=case_data.functional_group,
        atoms=case_data.atoms,
    )


def _test_get_atoms(functional_group, atoms):
    fg_atoms = it.zip_longest(functional_group.get_atoms(), atoms)
    for atom1, atom2 in fg_atoms:
        atom1 is atom2
