from .utilities import is_equivalent_atom


def test_clone(atom):
    """
    Test :meth:`.Atom.clone`.

    Parameters
    ----------
    atom : :class:`.Atom`
        The atom to be cloned.

    Returns
    -------
    None : :class:`NoneType`

    """

    clone = atom.clone()
    is_equivalent_atom(atom, clone)
