import pytest

from .utilities import is_clone


@pytest.mark.skip
def test_clone(construction_state):
    clone = construction_state.clone()
    is_clone(construction_state, clone)
