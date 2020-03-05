from .utilities import is_clone


def test_clone(case_data):
    vertex = case_data.vertex
    clone = vertex.clone()
    is_clone(vertex, clone)
