def check_bonds(old_state, new_state, building_blocks):
    old_bonds = tuple(old_state.get_bonds())
    new_bonds = tuple(new_state.get_bonds())
    num_added = len(new_bonds) - len(old_bonds)
    assert num_added == sum(bb.get_num_bonds() for bb in building_blocks)
    for bond1, bond2 in zip(
        new_bonds[len(old_bonds) :],
        get_expected_bonds(old_state, building_blocks),
    ):
        is_equivalent_bond(bond1, bond2)


def get_expected_bonds(old_state, building_blocks):
    start_atom_id = sum(1 for _ in old_state.get_atoms())
    for building_block in building_blocks:
        atom_map = {
            atom.get_id(): atom.with_id(id_)
            for id_, atom in enumerate(
                building_block.get_atoms(),
                start_atom_id,
            )
        }
        start_atom_id += building_block.get_num_atoms()
        for bond in building_block.get_bonds():
            yield bond.with_atoms(atom_map)


def is_equivalent_bond(bond1, bond2):
    is_equivalent_atom(bond1.get_atom1(), bond2.get_atom1())
    is_equivalent_atom(bond1.get_atom2(), bond2.get_atom2())


def is_equivalent_atom(atom1, atom2):
    assert atom1.__class__ is atom2.__class__
    assert atom1.get_id() == atom2.get_id()
    assert atom1.get_charge() == atom2.get_charge()
