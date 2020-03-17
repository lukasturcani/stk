from ..utilities import is_clone_functional_group


def test_clone(functional_group):
    """
    Test :meth:`.FunctionalGroup.clone`.

    Parameters
    ----------
    functional_group : :class:`.FunctionalGroup`
        The functional group to clone.

    Returns
    -------
    None : :class:`NoneType`

    """

    clone = functional_group.clone()
    is_clone_functional_group(functional_group, clone)
