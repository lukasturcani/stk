from stk import AtomInfo


def check_atom_infos(old_state, new_state, building_blocks):
    old_infos = tuple(old_state.get_atom_infos())
    new_infos = tuple(new_state.get_atom_infos())
    num_added = len(new_infos) - len(old_infos)
    assert num_added == sum(bb.get_num_atoms() for bb in building_blocks)

    for info1, info2 in zip(
        # Take just the newly added infos.
        new_infos[len(old_infos) :],
        get_infos(
            building_blocks=building_blocks,
            building_block_index=max(
                (info.building_block_index + 1 for info in old_infos),
                default=0,
            ),
            atom_id=len(old_infos),
        ),
    ):
        is_equivalent_atom(info1.atom, info2.atom)
        assert info1.building_block is info2.building_block
        assert info1.building_block_index is info2.building_block_index


def get_infos(building_blocks, building_block_index, atom_id):
    for building_block in building_blocks:
        for atom in building_block.get_atoms():
            yield AtomInfo(
                atom=atom.with_id(atom_id),
                building_block=building_block,
                building_block_index=building_block_index,
            )
            atom_id += 1
        building_block_index += 1


def is_equivalent_atom(atom1, atom2):
    assert atom1.__class__ is atom2.__class__
    assert atom1.get_id() == atom2.get_id()
    assert atom1.get_charge() == atom2.get_charge()
