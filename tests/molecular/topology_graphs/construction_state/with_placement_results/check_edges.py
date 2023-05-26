def check_edges(old_state, new_state):
    assert old_state.get_num_edges() == new_state.get_num_edges()
    for edge_id in range(old_state.get_num_edges()):
        assert old_state.get_edge(edge_id) is new_state.get_edge(edge_id)
