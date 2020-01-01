import itertools as it


def is_atom_clone(atom, clone):
    """
    Test that `clone` is a clone of `atom`.

    """

    assert atom is not clone
    assert atom.id == clone.id
    assert atom.charge == clone.charge
    assert atom.__class__ is clone.__class__


def is_clone_sequence(atoms, clones, atom_map):
    """
    Test if atoms are correctly cloned.

    """

    if atom_map is None:
        atom_map = {}

    for atom, clone in it.zip_longest(atoms, clones):
        is_atom_clone(atom_map.get(atom.id, atom), clone)


def is_clone_functional_group(
    functional_group1,
    functional_group2,
    atom_map,
):
    is_clone_sequence(
        functional_group1.get_atoms(),
        functional_group2.get_atoms(),
        atom_map,
    )
    is_clone_sequence(
        functional_group1.get_bonders(),
        functional_group2.get_bonders(),
        atom_map,
    )
    is_clone_sequence(
        functional_group1.get_deleters(),
        functional_group2.get_deleters(),
        atom_map,
    )
