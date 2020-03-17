import stk


def test_repr(bond):
    other = eval(repr(bond), dict(stk.__dict__))
    is_equivalent_atom(other.get_atom1(), bond.get_atom1())
    is_equivalent_atom(other.get_atom2(), bond.get_atom2())
    assert other.get_order() == bond.get_order()
    assert other.get_periodicity() == bond.get_periodicity()


def is_equivalent_atom(atom1, atom2):
    assert atom1.get_id() == atom2.get_id()
    assert atom1.get_charge() == atom2.get_charge()
    assert atom1.__class__ is atom2.__class__
    assert atom1.get_atomic_number() == atom2.get_atomic_number()
