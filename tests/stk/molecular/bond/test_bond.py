def test_clone(bond, atom_map):
    bond.attr = 1
    bond._attr = 2
    clone = bond.clone()
    assert bond.atom1 is not clone.atom1
    assert bond.atom1.id == clone.atom1.id
    assert bond.atom1.__class__ is clone.atom1.__class__
    assert bond.atom2 is not clone.atom2
    assert bond.atom2.id == clone.atom2.id
    assert bond.atom1.__class__ is clone.atom2.__class__
    assert bond.periodicity == clone.periodicity
    assert bond.attr == clone.attr
    assert not hasattr(clone, '_attr')


def test_is_periodic(bond):
    is_periodic = any(x != 0 for x in bond.periodicity)
    assert is_periodic == bond.is_periodic()
