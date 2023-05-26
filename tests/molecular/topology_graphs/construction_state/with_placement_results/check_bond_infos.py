from stk import BondInfo


def check_bond_infos(old_state, new_state, building_blocks):
    old_infos = tuple(old_state.get_bond_infos())
    new_infos = tuple(new_state.get_bond_infos())
    num_added = len(new_infos) - len(old_infos)
    assert num_added == sum(bb.get_num_bonds() for bb in building_blocks)

    for info1, info2 in zip(
        # Take just the newly added infos.
        new_infos[len(old_infos) :],
        get_infos(old_state, new_state, building_blocks),
    ):
        is_equivalent_bond(info1.bond, info2.bond)
        assert info1.building_block is info2.building_block
        assert info1.building_block_index is info2.building_block_index


def get_infos(old_state, new_state, building_blocks):
    start_building_block_index = max(
        (info.building_block_index + 1 for info in old_state.get_bond_infos()),
        default=0,
    )
    start_atom_index = sum(1 for _ in old_state.get_atoms())
    for building_block_index, building_block in enumerate(
        building_blocks,
        start_building_block_index,
    ):
        atom_map = {
            atom.get_id(): atom.with_id(id_)
            for id_, atom in enumerate(
                building_block.get_atoms(),
                start_atom_index,
            )
        }
        start_atom_index += building_block.get_num_atoms()

        for bond in building_block.get_bonds():
            yield BondInfo(
                bond=bond.with_atoms(atom_map),
                building_block=building_block,
                building_block_index=building_block_index,
            )


def is_equivalent_bond(bond1, bond2):
    is_equivalent_atom(bond1.get_atom1(), bond2.get_atom1())
    is_equivalent_atom(bond1.get_atom2(), bond2.get_atom2())


def is_equivalent_atom(atom1, atom2):
    assert atom1.__class__ is atom2.__class__
    assert atom1.get_id() == atom2.get_id()
    assert atom1.get_charge() == atom2.get_charge()
