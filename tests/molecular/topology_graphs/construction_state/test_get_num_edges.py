def test_get_num_edges(test_case):
    _test_get_num_edges(
        construction_state=test_case.construction_state,
        num_edges=test_case.get_num_edges(),
    )


def _test_get_num_edges(construction_state, num_edges):
    assert construction_state.get_num_edges() == num_edges
