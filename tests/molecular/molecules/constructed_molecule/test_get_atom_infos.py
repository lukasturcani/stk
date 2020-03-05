def test_get_atom_infos(case_data):
    _test_get_atom_infos(
        constructed_molecule=case_data.constructed_molecule,
        num_new_atoms=case_data.num_new_atoms,
    )


def _test_get_atom_infos(constructed_molecule, num_new_atoms):
    new_atoms = filter(
        lambda atom_info: atom_info.building_block is None,
        constructed_molecule.get_atom_infos(),
    )
    assert sum(1 for _ in new_atoms) == num_new_atoms
    assert (
        constructed_molecule.get_num_atoms()
        == sum(1 for _ in constructed_molecule.get_atom_infos())
    )
