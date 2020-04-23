import stk

from ..utilities import is_equivalent_functional_group


def test_repr(functional_group):
    """
    Test :meth:`.FunctionalGroup.__repr__`.

    Parameters
    ----------
    functional_group : :class:`.FunctionalGroup`
        The functional group whose representation is tested.

    Returns
    -------
    None : :class:`NoneType`

    """

    other = eval(repr(functional_group), dict(stk.__dict__))
    is_equivalent_functional_group(functional_group, other)
