import numpy as np


def is_clone(edge1, edge2):
    assert edge1.get_vertex1_id() == edge2.get_vertex1_id()
    assert edge1.get_vertex2_id() == edge2.get_vertex2_id()
    assert np.all(np.equal(edge1.get_position(), edge2.get_position()))
    assert edge1.get_periodicity() == edge2.get_periodicity()
    assert edge1.get_id() == edge2.get_id()
