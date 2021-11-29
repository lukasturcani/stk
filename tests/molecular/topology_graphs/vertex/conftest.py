import numpy as np
import pytest


class ScaleData:
    """
    Parameters for a :meth:`.Vertex.with_scale` test.

    Attributes
    ----------
    start : :class:`numpy.ndarray`
        The starting position of the tested vertex.

    scale : :class:`float` or :class:`numpy.ndarray`
        The scale applied to tested vertex.

    target : :class:`numpy.ndarray`
        The position of the new vertex.

    """

    def __init__(self, start, scale, target):
        """
        Initialize a :class:`.ScaleData` instance.

        Parameters
        ----------
        start : :class:`numpy.ndarray`
            The starting position of the tested vertex.

        scale : :class:`float` or :class:`numpy.ndarray`
            The scale applied to tested vertex.

        target : :class:`numpy.ndarray`
            The position of the new vertex.

        """

        self.start = start
        self.scale = scale
        self.target = target


@pytest.fixture(
    params=(
        ScaleData(np.array([0, 0, 0]), 2, np.array([0, 0, 0])),
        ScaleData(np.array([0, 0, 0]), 0, np.array([0, 0, 0])),
        ScaleData(np.array([1, -1, 0]), 0, np.array([0, 0, 0])),
        ScaleData(np.array([1, -1, 2]), 2, np.array([2, -2, 4])),
        ScaleData(np.array([2, 4, -5]), -3, np.array([-6, -12, 15])),
        ScaleData(
            np.array([1, -1, 2]),
            np.array([10, 2, -2]),
            np.array([10, -2, -4]),
        ),
    ),
)
def scale_data(request):
    """
    A :class:`.ScaleData` instance.

    """

    return request.param


@pytest.fixture(
    params=(
        [0, 0, 0],
        [1, 2, -20],
    ),
)
def position(request):
    return np.array(request.param, dtype=np.float64)


@pytest.fixture
def vertex(case_data):
    return case_data.vertex
