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
    assert clone.get_periodicity() is bond.get_periodicity()
    assert clone.get_order() is bond.get_order()
