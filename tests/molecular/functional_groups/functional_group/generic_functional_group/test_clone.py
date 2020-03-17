from ..utilities import is_clone_generic_functional_group


def test_clone(generic_functional_group):
    """
    Test :meth:`.GenericFunctionalGroup.clone`.

    Parameters
    ----------
    generic_functional_group : :class:`.GenericFunctionalGroup`
        The functional group to test.

    Returns
    -------
    None : :class:`NoneType`

    """

    clone = generic_functional_group.clone()
    is_clone_generic_functional_group(generic_functional_group, clone)
