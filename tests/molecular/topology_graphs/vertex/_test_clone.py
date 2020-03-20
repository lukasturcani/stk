from .utilities import is_clone


def test_clone(vertex):
    """
    Test :meth:`.Vertex.clone`.

    Parameters
    ----------
    vertex : :class:`.Vertex`
        The vertex to test.

    Returns
    -------
    None : :class:`NoneType`


    """

    clone = vertex.clone()
    is_clone(vertex, clone)
