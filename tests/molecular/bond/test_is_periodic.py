import stk


def test_is_periodic(atom1, atom2, order, periodicity):
    bond = stk.Bond(atom1, atom2, order, periodicity)
    return bond.is_periodic() == any(p != 0 for p in periodicity)
