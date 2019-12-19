import stk


def test_clone(atom):
    atom.attr = 1
    atom._attr = 2
    clone = atom.clone()
    assert clone is not atom
    assert clone.__class__ is atom.__class__
    assert clone.id == atom.id
    assert clone.charge == atom.charge
    assert clone.attr == atom.attr
    assert not hasattr(clone, '_attr')


def test_str(atom):
    atom.attr = 1
    atom._attr = 2
    other = eval(str(atom), dict(stk.__dict__))
    assert other is not atom
    assert other.__class__ is atom.__class__
    assert other.id == atom.id
    assert other.charge == atom.charge
    assert not hasattr(other, 'attr')
    assert not hasattr(other, '_attr')


def test_repr(atom):
    atom.attr = 1
    atom._attr = 2
    other = eval(repr(atom), dict(stk.__dict__))
    assert other is not atom
    assert other.__class__ is atom.__class__
    assert other.id == atom.id
    assert other.charge == atom.charge
    assert other.attr == atom.attr
    assert not hasattr(other, '_attr')
