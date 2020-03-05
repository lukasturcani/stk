import itertools as it


def test_get_atom_ids(case_data):
    _test_get_atom_ids(
        functional_group=case_data.functional_group,
        atoms=case_data.atoms,
    )


def _test_get_atom_ids(functional_group, atoms):
    fg_atoms = it.zip_longest(functional_group.get_atom_ids(), atoms)
    for id_, atom in fg_atoms:
        assert id_ == atom.get_id()
