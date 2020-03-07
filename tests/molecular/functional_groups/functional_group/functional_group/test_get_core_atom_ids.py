import itertools as it


def test_get_core_atom_ids(case_data):
    _test_get_core_atom_ids(
        functional_group=case_data.functional_group,
        core_atoms=case_data.core_atoms,
    )


def _test_get_core_atom_ids(functional_group, core_atoms):
    core_atoms = it.zip_longest(
        functional_group.get_core_atom_ids(),
        core_atoms,
    )
    for id_, core_atom in core_atoms:
        assert id_ == core_atom.get_id()
