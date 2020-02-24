from .utilities import is_clone


def test_clone(test_case):
    edge = test_case.edge
    clone = edge.clone()
    is_clone(edge, clone)
