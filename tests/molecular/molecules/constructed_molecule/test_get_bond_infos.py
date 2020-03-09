def test_get_bond_infos(case_data):
    _test_get_bond_infos(
        constructed_molecule=case_data.constructed_molecule,
        num_new_bonds=case_data.num_new_bonds,
    )


def _test_get_bond_infos(constructed_molecule, num_new_bonds):
    new_bonds = filter(
        lambda bond_info: bond_info.get_building_block() is None,
        constructed_molecule.get_bond_infos(),
    )
    assert sum(1 for _ in new_bonds) == num_new_bonds
    assert (
        constructed_molecule.get_num_bonds()
        == sum(1 for _ in constructed_molecule.get_bond_infos())
    )
