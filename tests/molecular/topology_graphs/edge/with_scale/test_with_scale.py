import numpy as np
import stk

from ..utilities import is_clone


def test_with_scale(case_data):
    """
    Test :meth:`.Edge.with_scale`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the starting position, the scale to apply
        and the target position.

    Returns
    -------
    None : :class:`NoneType`

    """

    edge = stk.Edge(
        id=0,
        vertex1=stk.Vertex(0, [0, 0, 0]),
        vertex2=stk.Vertex(1, [10, 0, 0]),
        position=case_data.start,
    )
    # Save a clone to check immutability.
    clone = edge.clone()
    _test_with_scale(edge, case_data.scale, case_data.target)
    is_clone(edge, clone)


def _test_with_scale(edge, scale, target):
    """
    Teste :meth:`.Edge.with_scale`.

    Parameters
    ----------
    edge : :class:`.Edge`
        The edge to test.

    scale : :class:`float` or :class:`numpy.ndarray`
        The scale to apply to `edge`.

    target : :class:`numpy.ndarray`
        The correct position of the new edge.

    Returns
    -------
    None : :class:`NoneType`

    """

    clone = edge.with_scale(scale)
    assert np.allclose(
        a=clone.get_position(),
        b=target,
        atol=1e-16,
    )
