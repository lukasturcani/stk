import stk

from .utilities import is_equivalent_atom


def test_repr(bond):
    """
    Test :meth:`.Bond.__repr__`.

    Parameters
    ----------
    bond : :class:`.Bond`
        The bond whose representation is tested.

    Returns
    -------
    None : :class:`NoneType`

    """

    other = eval(repr(bond), dict(stk.__dict__))
    is_equivalent_atom(other.get_atom1(), bond.get_atom1())
    is_equivalent_atom(other.get_atom2(), bond.get_atom2())
    assert other.get_order() == bond.get_order()
    assert other.get_periodicity() == bond.get_periodicity()
