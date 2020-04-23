import numpy as np

from .utilities import is_clone


def test_with_scale(vertex, scale_data):
    """
    Test :meth:`.Vertex.with_scale`.

    Parameters
    ----------
    vertex : :class:`.Vertex`
        The vertex to test.

    scale_data : :class:`.ScaleData`
        The parameter for :meth:`.Vertex.with_scale` to test.

    Returns
    -------
    None : :class:`NoneType`

    """

    # Save a clone to check immutability.
    clone = vertex.clone()
    _test_with_scale(
        vertex=vertex,
        start=scale_data.start,
        scale=scale_data.scale,
        target=scale_data.target,
    )
    is_clone(vertex, clone)


def _test_with_scale(vertex, start, scale, target):
    """
    Test :meth:`.Vertex.with_scale`.

    Parameters
    ----------
    vertex : :class:`.Vertex`
        The vertex to test.

    start : :class:`numpy.ndarray`
        The starting position of `vertex`.

    scale : :class:`float` or :class:`numpy.ndarray`
        The scale to apply to `vertex`.

    target : :class:`numpy.ndarray`
        The correct position of the new vertex.

    Returns
    -------
    None : :class:`NoneType`

    """

    clone = vertex.with_position(start).with_scale(scale)
    assert np.allclose(
        a=clone.get_position(),
        b=target,
        atol=1e-16,
    )
