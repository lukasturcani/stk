import stk
import numpy as np


def _get_periodic_position(cof, vertex1, vertex2, edge):
    assert vertex1 in edge.vertices and vertex2 in edge.vertices
    if vertex1 is vertex2:
        return vertex1.get_position()

    direction = 1 if vertex1 is edge.vertices[0] else -1
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
            periodic = any(dim != 0 for dim in edge.get_periodicity())
            edge_position = edge.get_position(vertex)
            if periodic:
                periodic_positions = (
                    _get_periodic_position(hexagonal, vertex, v2, edge)
                    for v2 in edge.vertices
                )
                expected = np.divide(
                    sum(periodic_positions),
                    len(edge.vertices)
                )
            else:
                expected = np.divide(
                    sum(v.get_position() for v in edge.vertices),
                    len(edge.vertices)
                )
            assert np.allclose(expected, edge_position, atol=1e-8)
