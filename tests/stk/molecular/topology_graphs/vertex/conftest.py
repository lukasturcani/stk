import numpy as np
import pytest
import stk


@pytest.fixture(
    params=,
)
def make_vertex(request):
    return request.param


@pytest.fixture(
    params=[

    ],
)
def make_edges(request):
    return request.param


@pytest.fixture(
    params=[0, 1, 0.5, 2],
)
def scale(request):
    return request.param


@pytest.fixture(
    params=[
        np.array([0, 0, 0]),
        np.array([10, 3, 5]),
        np.array([-1, -2, 0]),
        np.array([-1, -20, -3]),
        np.array([-23, 34, 12]),
    ],
)
def position(request):
    return np.array(request.param)


@pytest.fixture(
    params=[

    ],
)
def edge_ids(request):
    return list(request.param)


@pytest.fixture(
    params=[
        (0, 0, 0),
        (1, 0, 0),
        (-1, 0, 0),
        (10, 3, -1),
    ],
)
def cell(request):
    return tuple(request.param)


class Graph:
    def __init__(
        self,
        vertex,
        vertices,
        edges,
        edge_centroid,
        edge_plane_normal,
    ):
        self.vertex = vertex
        self.vertices = vertices
        self.edges = edges
        self.edge_centroid = edge_centroid
        self.edge_plane_normal = edge_plane_normal


@pytest.fixture(
    params=[

    ],
)
def graph(request):
    return request.param
