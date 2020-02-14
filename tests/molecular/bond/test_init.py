import stk


def test_init(atom1, atom2, order, periodicity):
    bond = stk.Bond(atom1, atom2, order, periodicity)
    assert bond.get_atom1() is atom1
    assert bond.get_atom2() is atom2
    assert bond.get_order() == order
    assert bond.get_periodicity() == periodicity
