from .utilities import is_clone


def test_clone(test_case):
    vertex = test_case.vertex
    clone = vertex.clone()
    is_clone(vertex, clone)
