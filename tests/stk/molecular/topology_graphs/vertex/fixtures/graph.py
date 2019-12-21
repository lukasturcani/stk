import stk
import numpy as np


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
    for i, edge in enumerate(edge_data):
        edge.id = i
    vertex = vertex_data.get_vertex()
    return _Graph(
        vertex=vertex,
        vertices=(vertex, ),
        edges=tuple(e.get_edge() for e in edge_data),
        edge_centroid=np.array([0., 0., 0.]),
        reference=np.array([0., 0., 1.]),
        edge_plane_normal=np.array([0., 0., 1.]),
    )


def graph2(make_vertex_data):
    vertex_data = make_vertex_data(0, 0, 0)
    edge_data = (
        stk.EdgeData(
            vertex_data,
            vertex_data,
            periodicity=(1, 0, 0),
            lattice_constants=(
                np.array([1., 0., 0.]),
                np.array([0., 1., 0.]),
                np.array([0., 0., 1.]),
            ),
        ),
        stk.EdgeData(
            vertex_data,
            vertex_data,
            periodicity=(-1, 0, 0),
            lattice_constants=(
                np.array([1., 0., 0.]),
                np.array([0., 1., 0.]),
                np.array([0., 0., 1.]),
            ),
        ),
    )
    for i, edge in enumerate(edge_data):
        edge.id = i
    vertex = vertex_data.get_vertex()
    return _Graph(
        vertex=vertex,
        vertices=(vertex, ),
        edges=tuple(e.get_edge() for e in edge_data),
        edge_centroid=np.array([0., 0., 0.]),
        reference=np.array([0., 0., 1.]),
        edge_plane_normal=np.array([0., 0., 1.]),
    )


def _graph3(vertex_data):
    pass
