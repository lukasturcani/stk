import stk

from .utilities import is_equivalent_atom


def test_repr(atom):
    """
    Test :meth:`.Atom.__repr__`.

    Parameters
    ----------
    atom : :class:`.Atom`
        The atom, whose representation should be tested.

    Returns
    -------
    None : :class:`NoneType`

    """

    other = eval(repr(atom), dict(stk.__dict__))
    is_equivalent_atom(other, atom)
