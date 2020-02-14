import itertools as it


def is_equivalent_atom(atom1, atom2):
    assert atom1.get_id() == atom2.get_id()
    assert atom1.get_charge() == atom2.get_charge()
    assert atom1.get_mass() == atom2.get_mass()
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
    is_clone_sequence(
        functional_group1.get_bonders(),
        functional_group2.get_bonders(),
    )
    is_clone_sequence(
        functional_group1.get_deleters(),
        functional_group2.get_deleters(),
    )


def is_equivalent_sequence(atoms1, atoms2):
    for atom1, atom2 in it.zip_longest(atoms1, atoms2):
        is_equivalent_atom(atom1, atom2)


def is_equivalent_functional_group(
    functional_group1,
    functional_group2,
):
    is_equivalent_sequence(
        functional_group1.get_atoms(),
        functional_group2.get_atoms(),
    )
    is_equivalent_sequence(
        functional_group1.get_bonders(),
        functional_group2.get_bonders(),
    )
    is_equivalent_sequence(
        functional_group1.get_deleters(),
        functional_group2.get_deleters(),
    )
