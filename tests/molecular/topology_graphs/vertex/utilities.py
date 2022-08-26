import numpy as np


def is_clone(vertex, clone):
    assert np.all(
        np.equal(
            vertex.get_position(),
            clone.get_position(),
        )
    )
    assert vertex.get_id() == clone.get_id()
    assert np.all(np.equal(vertex.get_cell(), clone.get_cell()))
