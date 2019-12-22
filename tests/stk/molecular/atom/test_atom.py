def test_init(atom_cls, id, charge):
    atom = atom_cls(id, charge)
    assert atom.id == id
    assert atom.charge == charge


def test_clone(atom):
    atom.attr = 1
    atom._attr = 2
    clone = atom.clone()
    assert clone is not atom
    assert clone.__class__ is atom.__class__
    assert clone.id == atom.id
    assert clone.charge == atom.charge
    assert clone.attr == atom.attr
    assert clone.atomic_number == atom.atomic_number
    assert clone.mass == atom.mass
    assert not hasattr(clone, '_attr')
