import stk
import numpy as np
import pytest


class _Graph:
    def __init__(
        self,
        vertex,
        vertices,
        edges,
        edge_centroid,
        reference,
        edge_plane_normal,
    ):
        self.vertex = vertex
        self.vertices = vertices
        self.edges = edges
        self.edge_centroid = edge_centroid
        self.reference = reference
        self.edge_plane_normal = edge_plane_normal


def set_ids(data):
    for i in enumerate(data):
        data.id = i


@pytest.fixture
def graph1(make_vertex_data):
    vertex_data = make_vertex_data(0, 0, 0)
    edge_data = (
        stk.EdgeData(
            vertex_data,
            position=np.array([1., 0., 0.]),
        ),
        stk.EdgeData(
            vertex_data,
            position=np.array([-1., 0., 0.]),
        ),
        stk.EdgeData(
            vertex_data,
            position=np.array([0., 1., 0.]),
        ),
        stk.EdgeData(
            vertex_data,
            position=np.array([0., -1., 0.]),
        ),
    )
    set_ids(edge_data)
    vertex = vertex_data.get_vertex()
    return _Graph(
        vertex=vertex,
        vertices=(vertex, ),
        edges=tuple(e.get_edge() for e in edge_data),
        edge_centroid=np.array([0., 0., 0.]),
        reference=np.array([0., 0., 1.]),
        edge_plane_normal=np.array([0., 0., 1.]),
    )


@pytest.fixture(
    params=[
        (0, np.array([0., 0., 0.])),
        (1, np.array([0., 0., 0.])),
    ]
)
def graph2(request, make_vertex_data):
    vertex_data = (
        make_vertex_data(0, 0, 0),
        make_vertex_data(1, 0, 0),
    )
    edge_data = (
        stk.EdgeData(
            vertex_data[0],
            vertex_data[1],
            periodicity=(1, 0, 0),
            lattice_constants=(
                np.array([1., 0., 0.]),
                np.array([0., 1., 0.]),
                np.array([0., 0., 1.]),
            ),
        ),
        stk.EdgeData(
            vertex_data[0],
            vertex_data[1],
        ),
    )
    set_ids(vertex_data)
    set_ids(edge_data)
    vertices = tuple(v.get_vertex() for v in vertex_data)

    vertex_id, edge_centroid = request.param
    return _Graph(
        vertex=vertices[vertex_id],
        vertices=vertices,
        edges=tuple(e.get_edge() for e in edge_data),
        edge_centroid=edge_centroid,
        reference=np.array([0., 0., 1.]),
        edge_plane_normal=np.array([0., 0., 1.]),
    )


def _graph3(vertex_data):
    pass
