def check_vertices(old_state, new_state):
    assert old_state.get_num_vertices() == new_state.get_num_vertices()
    for vertex_id in range(old_state.get_num_vertices()):
        assert old_state.get_vertex(vertex_id) is new_state.get_vertex(
            vertex_id
        )
