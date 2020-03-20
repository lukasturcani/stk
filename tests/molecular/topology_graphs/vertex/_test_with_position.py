import numpy as np

from .utilities import is_clone


def test_with_position(vertex, position):
    """
    Test :meth:`.Vertex.with_position`.

    Parameters
    ----------
    vertex : :class:`.Vertex`
        The vertex to test.

    position : :class:`numpy.ndarray`
        The position of the new vertex.

    Returns
    -------
    None : :class:`NoneType`

    """

    # Save a clone to check immutability.
    clone = vertex.clone()
    _test_with_position(vertex, position)
    is_clone(vertex, clone)


def _test_with_position(vertex, position):
    """
    Test :meth:`.Vertex.with_position`.

    Parameters
    ----------
    vertex : :class:`.Vertex`
        The vertex to test.

    position : :class:`numpy.ndarray`
        The position of the new vertex.

    Returns
    -------
    None : :class:`NoneType`

    """

    clone = vertex.with_position(position)
    assert np.all(np.equal(clone.get_position(), position))
