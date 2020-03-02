from .utilities import is_clone


def test_clone(construction_state):
    clone = construction_state.clone()
    is_clone(construction_state, clone)
