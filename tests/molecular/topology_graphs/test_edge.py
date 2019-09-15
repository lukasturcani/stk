import stk
import numpy as np


def _get_periodic_position(cof, vertex1, vertex2, edge):
    edge_vertices = {
        cof.vertices[v_id] for v_id in edge.get_vertex_ids()
    }
    assert vertex1 in edge_vertices and vertex2 in edge_vertices
    if vertex1 is vertex2:
        return vertex1.get_position()

    v1_id, v2_id = list(edge.get_vertex_ids())
    direction = 1 if vertex1 is cof.vertices[v1_id] else -1
    cell_end = vertex1.get_cell() + direction*edge.get_periodicity()
    cell_shift = cell_end - vertex2.get_cell()
    shift = 0
    for dim, constant in zip(cell_shift, cof._lattice_constants):
        shift += dim*constant
    return vertex2.get_position() + shift


def test_get_position():
    hexagonal = stk.cof.Hexagonal((2, 2, 1))
    for vertex in hexagonal.vertices:
        for edge_id in vertex.get_edge_ids():
            edge = hexagonal.edges[edge_id]
            edge_position = edge.get_position(
                reference=vertex,
                vertices=hexagonal.vertices
            )
            if edge.is_periodic():
                periodic_positions = (
                    _get_periodic_position(
                        cof=hexagonal,
                        vertex1=vertex,
                        vertex2=hexagonal.vertices[v2_id],
                        edge=edge
                    )
                    for v2_id in edge.get_vertex_ids()
                )
                expected = np.divide(
                    sum(periodic_positions),
                    sum(1 for _ in edge.get_vertex_ids())
                )
            else:
                expected = np.divide(
                    sum(
                        hexagonal.vertices[v].get_position()
                        for v in edge.get_vertex_ids()
                    ),
                    sum(1 for _ in edge.get_vertex_ids())
                )
            assert np.allclose(expected, edge_position, atol=1e-8)
