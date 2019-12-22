import numpy as np


def test_apply_scale(vertex_data, position, scale):
    vertex_data.position = position
    vertex = vertex_data.get_vertex()
    before = vertex.get_position()
    vertex.apply_scale(scale)
    assert np.allclose(
        a=vertex.get_position(),
        b=scale*before,
        atol=1e-32,
    )


def test_clone(vertex_data):
    vertex = vertex_data.get_vertex()
    vertex.attr = 1
    vertex._attr = 2
    clone = vertex.clone()
    assert (
        np.all(np.equal(clone.get_position(), vertex.get_position()))
    )
    assert clone.attr == vertex.attr
    assert not hasattr(clone, '_attr')
    assert clone.get_num_edges() == vertex.get_num_edges()
    assert list(clone.get_edge_ids()) == list(vertex.get_edge_ids())
    assert np.all(np.equal(clone.get_cell(), vertex.get_cell()))


def test_get_position(vertex_data, position):
    vertex_data.position = position
    position = vertex_data.get_vertex().get_position()
    assert np.all(np.equal(position, vertex_data.position))


def test_get_num_edges(vertex_data, add_edge_data):
    add_edge_data(vertex_data)
    vertex = vertex_data.get_vertex()
    assert vertex.get_num_edges() == len(vertex_data.edges)


def test_get_edge_ids(vertex_data, add_edge_data):
    add_edge_data(vertex_data)
    vertex = vertex_data.get_vertex()
    assert (
        tuple(vertex.get_edge_ids())
        == tuple(data.id for data in vertex_data.edges)
    )


def test_get_cell(vertex_data):
    vertex = vertex_data.get_vertex()
    assert np.all(np.equal(vertex.get_cell(), vertex_data.cell))


def test_set_constructed_molecule(
    vertex_data,
    constructed_molecule,
    get_atom_ids,
):
    vertex = vertex_data.get_vertex()
    vertex.set_constructed_molecule(constructed_molecule)
    atom_ids = get_atom_ids(constructed_molecule)
    centroid = constructed_molecule.get_centroid(atom_ids)
    constructed_molecule._position_matrix = (
        constructed_molecule._position_matrix.T
    )
    assert np.allclose(
        a=vertex._get_molecule_centroid(atom_ids=atom_ids),
        b=centroid,
        atol=1e-15,
    )


def test_get_edge_centroid(centroid_graph):
    centroid_edges = tuple(
        centroid_graph.edges[edge_id]
        for edge_id in centroid_graph.vertex.get_edge_ids()
    )
    centroid = centroid_graph.vertex._get_edge_centroid(
        centroid_edges=centroid_edges,
        vertices=centroid_graph.vertices,
    )
    assert np.allclose(
        a=centroid,
        b=centroid_graph.edge_centroid,
        atol=1e-32,
    )


def test_get_edge_plane_normal(plane_graph):
    plane_edges = tuple(
        plane_graph.edges[edge_id]
        for edge_id in plane_graph.vertex.get_edge_ids()
    )
    normal = plane_graph.vertex._get_edge_plane_normal(
        reference=plane_graph.reference,
        plane_edges=plane_edges,
        vertices=plane_graph.vertices,
    )
    assert np.allclose(
        a=normal,
        b=plane_graph.edge_plane_normal,
        atol=1e-32,
    )


def test_get_molecule_centroid(
    vertex_data,
    constructed_molecule,
    get_atom_ids,
):
    vertex = vertex_data.get_vertex()
    atom_ids = get_atom_ids(constructed_molecule)
    constructed_molecule.set_centroid(vertex_data.position, atom_ids)
    constructed_molecule._position_matrix = (
        constructed_molecule._position_matrix.T
    )
    vertex.set_constructed_molecule(constructed_molecule)
    assert np.allclose(
        a=vertex._get_molecule_centroid(atom_ids),
        b=vertex_data.position,
        atol=1e-15,
    )
