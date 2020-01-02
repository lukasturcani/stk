import itertools as it


def test_get_atom_ids(test_case):
    _test_get_atom_ids(
        functional_group=test_case.functional_group,
        atoms=test_case.atoms,
    )


def _test_get_atom_ids(functional_group, atoms):
    fg_atoms = it.zip_longest(functional_group.get_atom_ids(), atoms)
    for id_, atom in fg_atoms:
        assert id_ == atom.get_id()
