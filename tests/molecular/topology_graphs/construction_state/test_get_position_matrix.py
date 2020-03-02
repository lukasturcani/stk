import numpy as np


def test_get_position_matrix(test_case):
    _test_get_position_matrix(
        construction_state=test_case.construction_state,
        position_matrix=test_case.position_matrix,
    )


def _test_get_position_matrix(construction_state, position_matrix):
    assert np.all(np.equal(
        construction_state.get_position_matrix(),
        position_matrix,
    ))
