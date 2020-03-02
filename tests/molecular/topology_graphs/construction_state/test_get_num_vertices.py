def test_get_num_vertices(test_case):
    _test_get_num_vertices(
        construction_state=test_case.construction_state,
        num_vertices=test_case.num_vertices,
    )


def _test_get_num_vertices(construction_state, num_vertices):
    assert construction_state.get_num_vertices() == num_vertices
