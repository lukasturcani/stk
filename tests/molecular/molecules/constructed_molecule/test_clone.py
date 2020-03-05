import itertools as it


def test_clone(case_data):
    _test_clone(case_data.constructed_molecule)


def _test_clone(constructed_molecule):
    clone = constructed_molecule.clone()
    is_clone(constructed_molecule, clone)


def is_clone(molecule1, molecule2):
    assert (
        molecule1.get_topology_graph()
        is molecule2.get_topology_graph()
    )
    assert (
        molecule1.get_building_block_counts()
        == molecule2.get_building_block_counts()
    )
    assert (
        molecule1.get_building_block_vertices()
        == molecule2.get_building_block_vertices()
    )
    for building_block1, building_block2 in it.zip_longest(
        molecule1.get_building_blocks(),
        molecule2.get_building_blocks(),
    ):
        assert building_block1 is building_block2
    for atom_info1, atom_info2 in it.zip_longest(
        molecule1.get_atom_infos(),
        molecule2.get_atom_infos(),
    ):
        assert atom_info1 is atom_info2
    for bond_info1, bond_info2 in it.zip_longest(
        molecule1.get_bond_infos(),
        molecule2.get_bond_infos(),
    ):
        assert bond_info1 is bond_info2
