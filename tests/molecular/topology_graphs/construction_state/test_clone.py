from .utilities import is_clone


def test_clone(test_case):
    construction_state = test_case.construction_state
    clone = construction_state.clone()
    is_clone(construction_state, clone)
