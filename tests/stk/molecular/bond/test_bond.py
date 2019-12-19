import stk


def test_clone(bond):
    bond.attr = 1
    bond._attr = 2
    clone = bond.clone()
    assert bond.atom1 is not clone.atom1
    assert bond.atom1.id == clone.atom1.id
    assert bond.atom2 is not clone.atom2
    assert bond.atom2.id == clone.atom2.id
    assert bond.periodicity == clone.periodicity
    assert bond.attr == clone.attr
    assert not hasattr(clone, '_attr')


def test_is_periodic(bond):
    is_periodic = any(x != 0 for x in bond.periodicity)
    assert is_periodic == bond.is_periodic()


def test_str(bond):
    bond.attr = 1
    bond._attr = 2
    other = eval(str(bond), dict(stk.__dict__))
    assert bond.atom1 is not other.atom1
    assert bond.atom1.id == other.atom1.id
    assert bond.atom2 is not other.atom2
    assert bond.atom2.id == other.atom2.id
    assert bond.periodicity == other.periodicity
    assert not hasattr(other, 'attr')
    assert not hasattr(other, '_attr')


def test_repr(bond):
    bond.attr = 1
    bond._attr = 2
    other = eval(repr(bond), dict(stk.__dict__))
    assert bond.atom1 is not other.atom1
    assert bond.atom1.id == other.atom1.id
    assert bond.atom2 is not other.atom2
    assert bond.atom2.id == other.atom2.id
    assert bond.periodicity == other.periodicity
    assert bond.attr == other.attr
    assert not hasattr(other, '_attr')
