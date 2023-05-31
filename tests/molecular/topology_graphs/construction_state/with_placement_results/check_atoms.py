import itertools as it


def check_atoms(old_state, new_state, building_blocks):
    old_atoms = tuple(old_state.get_atoms())
    new_atoms = tuple(new_state.get_atoms())
    num_added = len(new_atoms) - len(old_atoms)
    assert num_added == sum(bb.get_num_atoms() for bb in building_blocks)
    expected_atoms = enumerate(
        it.chain(*(bb.get_atoms() for bb in building_blocks)),
        len(old_atoms),
    )
    for atom1, (id_, atom2) in zip(
        # Take just the newly added atoms.
        new_atoms[len(old_atoms) :],
        expected_atoms,
    ):
        is_equivalent_atom(atom1, atom2.with_id(id_))


def is_equivalent_atom(atom1, atom2):
    assert atom1.get_id() == atom2.get_id()
    assert atom1.__class__ is atom2.__class__
    assert atom1.get_charge() == atom2.get_charge()
