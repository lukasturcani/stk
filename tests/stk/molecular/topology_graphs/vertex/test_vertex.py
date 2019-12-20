import numpy as np


def test_apply_scale(make_vertex, position, scale):
    vertex = make_vertex(position)
    before = vertex.get_position()
    vertex.apply_scale(scale)
    assert np.allclose(
        a=vertex.get_position(),
        b=scale*before,
        atol=1e-32,
    )


def test_clone(make_vertex):
    vertex = make_vertex()
    vertex.attr = 1
    vertex._attr = 2
    clone = vertex.clone()
    assert clone.get_position() == vertex.get_position()
    assert clone.attr == vertex.attr
    assert not hasattr(clone, '_attr')
    assert clone.get_num_edges() == vertex.get_num_edges()
    assert list(clone.get_edge_ids()) == list(vertex.get_edge_ids())
    assert clone.get_cell() == vertex.get_cell()


def test_get_position(make_vertex, position):
    vertex = make_vertex(position)
    assert np.all(np.equal(vertex.get_position(), position))


def test_get_num_edges(make_vertex, make_edges, edge_ids):
    vertex = make_vertex(edges=make_edges(edge_ids))
    assert vertex.get_num_edges() == len(edge_ids)


def test_get_edge_ids(make_vertex, make_edges, edge_ids):
    vertex = make_vertex(edges=make_edges(edge_ids))
    assert tuple(vertex.get_edge_ids()) == edge_ids


def test_get_cell(make_vertex, cell):
    vertex = make_vertex(cell=cell)
    assert vertex.get_cell() == cell


def test_set_constructed_molecule(make_vertex, constructed_molecule):
    vertex = make_vertex()
    vertex.set_constructed_molecule(molecule)
    assert


def test_get_edge_centroid(vertex_test_case, edges):
    vertex = vertex_test_case.vertex
    centroid_edges = vertex_test_case.centroid_edges
    vertices = vertex_test_case.vertices
    assert (
        vertex._get_edge_centroid(centroid_edges, vertices)
        == vertex_test_case.edge_centroid
    )


def test_get_edge_plane_normal(vertex_test_case, edges):
    vertex = vertex_test_case.vertex


def test_get_molecule_centroid(vertex_test_case, edges, molecule):
    assert False
