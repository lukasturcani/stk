def test_clone(bond):
    bond.attr = 1
    bond._attr = 2
    clone = bond.clone()

    assert clone.get_atom1() is bond.get_atom1()
    assert clone.get_atom2() is bond.get_atom2()
    assert bond.get_periodicity() == clone.get_periodicity()
    assert bond.get_order() == clone.get_order()
    assert bond.attr == clone.attr
    assert not hasattr(clone, '_attr')
