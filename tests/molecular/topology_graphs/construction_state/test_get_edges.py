import itertools as it


def test_get_edges(test_case):
    _test_get_edges(
        construction_state=test_case.construction_state,
        vertex_edges=test_case.vertex_edges,
    )


def _test_get_edges(construction_state, vertex_edges):
    for vertex_id, edge in vertex_edges.items():
        assert construction_state.get_vertex_edges(vertex_id) is edge
