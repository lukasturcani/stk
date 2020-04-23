from .utilities import is_clone


def test_clone(edge):
    """
    Test :meth:`.Edge.clone`.

    Parameters
    ----------
    edge : :class:`.Edge`
        The edge to test.

    Returns
    -------
    None : :class:`NoneType`

    """

    clone = edge.clone()
    is_clone(edge, clone)
