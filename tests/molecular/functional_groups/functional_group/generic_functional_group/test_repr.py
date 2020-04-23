import stk

from ..utilities import is_equivalent_generic_functional_group


def test_repr(generic_functional_group):
    """
    Test :meth:`.GenericFunctionalGroup.__repr__`.

    Parameters
    ----------
    generic_functional_group : :class:`.GenericFunctionalGroup`
        The functional group, whose representation should be tested.

    Returns
    -------
    None : :class:`NoneType`

    """

    other = eval(repr(generic_functional_group), dict(stk.__dict__))
    is_equivalent_generic_functional_group(
        functional_group1=generic_functional_group,
        functional_group2=other,
    )
