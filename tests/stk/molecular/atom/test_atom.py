def test_init(get_atom, id, charge):
    atom = get_atom(id, charge)
    assert atom.get_id() == id
    assert atom.get_charge() == charge


def test_clone(atom):
    atom.attr = 1
    atom._attr = 2
    clone = atom.clone()
    assert clone is not atom
    assert clone.__class__ is atom.__class__
    assert clone.get_id() == atom.get_id()
    assert clone.get_charge() == atom.get_charge()
    assert clone.attr == atom.attr
    assert clone.get_atomic_number() == atom.get_atomic_number()
    assert clone.get_mass() == atom.get_mass()
    assert not hasattr(clone, '_attr')
