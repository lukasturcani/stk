from .utilities import is_equivalent_atom


def test_with_id(atom, id):
    """
    Test :meth:`.Atom.with_id`.

    Parameters
    ----------
    atom : :class:`.Atom`
        The atom to test.

    id : :class:`int`
        The correct id of the new atom.

    Returns
    -------
    None : :class:`NoneType`

    """

    # Save a clone to ensure that "atom" is not changed by the test.
    before = atom.clone()
    _test_with_id(atom, id)
    is_equivalent_atom(before, atom)


def _test_with_id(atom, id):
    """
    Test :meth:`.Atom.with_id`.

    Parameters
    ----------
    atom : :class:`.Atom`
        The atom to test.

    id : :class:`int`
        The correct id of the new atom.

    Returns
    -------
    None : :class:`NoneType`

    """

    new_atom = atom.with_id(id)
    assert new_atom is not atom
    assert new_atom.__class__ is atom.__class__
    assert new_atom.get_id() == id
    assert new_atom.get_charge() == atom.get_charge()
    assert new_atom.get_atomic_number() == atom.get_atomic_number()
