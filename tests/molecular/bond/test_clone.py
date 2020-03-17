def test_clone(bond):
    """
    Test :meth:`.Bond.clone`.

    Parameters
    ----------
    bond : :class:`.Bond`
        The bond to test.

    Returns
    -------
    None : :class:`NoneType`

    """

    clone = bond.clone()

    assert clone.get_atom1() is bond.get_atom1()
    assert clone.get_atom2() is bond.get_atom2()
    assert bond.get_periodicity() is clone.get_periodicity()
    assert bond.get_order() is clone.get_order()
