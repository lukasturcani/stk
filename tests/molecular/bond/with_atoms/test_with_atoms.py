def test_with_atoms(bond, get_atom_map):
    """
    Test :meth:`.Bond.with_atoms`.

    Parameters
    ----------
    bond : :class:`.Bond`
        The bond to test.

    get_atom_map : :class:`callable`
        Takes a single parameter, `bond`, and returns a valid
        `atom_map` parameter for its :meth:`.Bond.with_atoms` method.
        This allows the testing of different values of this parameter.

    Returns
    -------
    None : :class:`NoneType`

    """

    atom_map = get_atom_map(bond)
    clone = bond.with_atoms(atom_map)
    assert clone is not bond

    expected_atom1 = atom_map.get(
        bond.get_atom1().get_id(),
        bond.get_atom1(),
    )
    assert clone.get_atom1() is expected_atom1

    expected_atom2 = atom_map.get(
        bond.get_atom2().get_id(),
        bond.get_atom2(),
    )
    assert clone.get_atom2() is expected_atom2

    assert bond.get_periodicity() == clone.get_periodicity()
    assert bond.get_order() == clone.get_order()
