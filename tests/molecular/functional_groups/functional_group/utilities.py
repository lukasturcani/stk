import itertools as it


def is_equivalent_atom(atom1, atom2):
    assert atom1.get_id() == atom2.get_id()
    assert atom1.get_charge() == atom2.get_charge()
    assert atom1.get_atomic_number() == atom2.get_atomic_number()
    assert atom1.__class__ is atom2.__class__


def is_clone_sequence(atoms, clones):
    """
    Test if atoms are correctly cloned.

    """

    for atom, clone in it.zip_longest(atoms, clones):
        assert atom is clone


def is_clone_functional_group(functional_group1, functional_group2):
    assert functional_group1 is not functional_group2

    is_clone_sequence(
        functional_group1.get_atoms(),
        functional_group2.get_atoms(),
    )
    same_ids(
        ids1=functional_group1.get_atom_ids(),
        ids2=functional_group2.get_atom_ids(),
    )
    same_ids(
        ids1=functional_group1.get_placer_ids(),
        ids2=functional_group2.get_placer_ids(),
    )
    same_ids(
        ids1=functional_group1.get_core_atom_ids(),
        ids2=functional_group2.get_core_atom_ids(),
    )


def is_clone_generic_functional_group(
    functional_group1,
    functional_group2,
):
    is_clone_functional_group(functional_group1, functional_group2)

    is_clone_sequence(
        functional_group1.get_bonders(),
        functional_group2.get_bonders(),
    )
    same_ids(
        functional_group1.get_bonder_ids(),
        functional_group2.get_bonder_ids(),
    )

    is_clone_sequence(
        functional_group1.get_deleters(),
        functional_group2.get_deleters(),
    )
    same_ids(
        functional_group1.get_deleter_ids(),
        functional_group2.get_deleter_ids(),
    )


def is_equivalent_sequence(atoms1, atoms2):
    for atom1, atom2 in it.zip_longest(atoms1, atoms2):
        is_equivalent_atom(atom1, atom2)


def same_ids(ids1, ids2):
    for id1, id2 in it.zip_longest(ids1, ids2):
        assert id1 == id2


def is_equivalent_functional_group(
    functional_group1,
    functional_group2,
):
    is_equivalent_sequence(
        functional_group1.get_atoms(),
        functional_group2.get_atoms(),
    )
    same_ids(
        ids1=functional_group1.get_atom_ids(),
        ids2=functional_group2.get_atom_ids(),
    )
    same_ids(
        ids1=functional_group1.get_placer_ids(),
        ids2=functional_group2.get_placer_ids(),
    )
    same_ids(
        ids1=functional_group1.get_core_atom_ids(),
        ids2=functional_group2.get_core_atom_ids(),
    )


def is_equivalent_generic_functional_group(
    functional_group1,
    functional_group2,
):
    is_equivalent_functional_group(
        functional_group1=functional_group1,
        functional_group2=functional_group2,
    )
    is_equivalent_sequence(
        functional_group1.get_bonders(),
        functional_group2.get_bonders(),
    )
    same_ids(
        ids1=functional_group1.get_bonder_ids(),
        ids2=functional_group2.get_bonder_ids(),
    )
    assert (
        functional_group1.get_num_bonders()
        == functional_group2.get_num_bonders()
    )

    is_equivalent_sequence(
        functional_group1.get_deleters(),
        functional_group2.get_deleters(),
    )
    same_ids(
        ids1=functional_group1.get_deleter_ids(),
        ids2=functional_group2.get_deleter_ids(),
    )
