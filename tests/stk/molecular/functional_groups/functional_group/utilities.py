def is_atom_clone(atom, clone):
    """
    Test that `clone` is a clone of `atom`.

    """

    assert atom is not clone
    assert atom.id == clone.id
    assert atom.charge == clone.charge
    assert atom.__class__ is clone.__class__
