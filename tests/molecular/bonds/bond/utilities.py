def is_equivalent_atom(atom1, atom2):
    assert atom1.get_id() == atom2.get_id()
    assert atom1.get_charge() == atom2.get_charge()
    assert atom1.__class__ is atom2.__class__
    assert atom1.get_atomic_number() == atom2.get_atomic_number()
