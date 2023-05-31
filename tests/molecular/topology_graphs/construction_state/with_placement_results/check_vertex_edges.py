def check_vertex_edges(old_state, new_state):
    for vertex_id in range(old_state.get_num_vertices()):
        assert old_state.get_edges(vertex_id) == new_state.get_edges(vertex_id)
