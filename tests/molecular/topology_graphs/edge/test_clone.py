from .utilities import is_clone


def test_clone(case_data):
    edge = case_data.edge
    clone = edge.clone()
    is_clone(edge, clone)
