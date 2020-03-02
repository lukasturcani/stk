import itertools as it
import numpy as np


def test_get_vertex(test_case):
    _test_get_vertex(
        construction_state=test_case.construction_state,
        vertices=test_case.vertices,
    )


def _test_get_vertex(construction_state, vertices):
    state_vertices = map(
        construction_state.get_vertex,
        range(construction_state.get_num_vertices()),
    )
    for vertex1, vertex2 in it.zip_longest(state_vertices, vertices):
        assert vertex1.get_id() == vertex2.get_id()
        assert np.all(np.equal(
            vertex1.get_position(),
            vertex2.get_position(),
        ))
        assert vertex1.__class__ is vertex2.__class__
