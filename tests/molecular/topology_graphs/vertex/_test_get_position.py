import numpy as np


def test_get_position(case_data):
    """
    Test :meth:`.Vertex.get_position`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the vertex to test and the correct
        position.

    Returns
    -------
    None : :class:`NoneType`


    """

    _test_get_position(case_data.vertex, case_data.position)


def _test_get_position(vertex, position):
    """
    Test :meth:`.Vertex.get_position`.

    Parameters
    ----------
    vertex : :class:`.Vertex`
        The vertex to test.

    position : :class:`numpy.ndarray`
        The correct position.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert np.allclose(vertex.get_position(), position, atol=1e-15)
