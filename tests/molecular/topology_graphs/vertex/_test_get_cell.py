import numpy as np


def test_get_cell(case_data):
    """
    Test :meth:`.Vertex.get_cell`.

    Parameters
    ----------
    case_data : :class:`.CaseData`
        A test case. Holds the vertex to test and the correct cell.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_get_cell(case_data.vertex, case_data.cell)


def _test_get_cell(vertex, cell):
    """
    Test :meth:`.Vertex.get_cell`.

    Parameters
    ----------
    vertex : :class:`.Vertex`
        The vertex to test.

    cell : :class:`numpy.ndarray`
        The correct cell.

    Returns
    -------
    None : :class:`NoneType`

    """

    assert np.all(np.equal(vertex.get_cell(), cell))
