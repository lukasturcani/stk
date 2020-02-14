def test_to_dict(bond):
    d = bond.to_dict()
    assert d['atom1_id'] == bond.get_atom1().get_id()
    assert d['atom2_id'] == bond.get_atom2().get_id()
    assert d['order'] == bond.get_order()
    assert tuple(d['periodicity']) == bond.get_periodicity()
