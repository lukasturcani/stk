import numpy as np


def check_position_matrix(old_state, new_state, placement_results):
    assert np.all(
        np.equal(
            get_expected_position_matrix(old_state, placement_results),
            new_state.get_position_matrix(),
        )
    )


def get_expected_position_matrix(old_state, placement_results):
    return np.vstack(
        [
            old_state.get_position_matrix(),
            *(result.position_matrix for result in placement_results),
        ]
    )
