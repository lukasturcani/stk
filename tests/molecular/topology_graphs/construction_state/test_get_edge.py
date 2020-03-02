import itertools as it
import numpy as np


def test_get_edge(test_case):
    _test_get_edge(
        construction_state=test_case.construction_state,
        edges=test_case.edges,
    )


def _test_get_edge(construction_state, edges):
    state_edges = map(
        construction_state.get_edge,
        range(construction_state.get_num_edges),
    )
    for edge1, edge2 in it.zip_longest(state_edges, edges):
        assert edge1.get_id() == edge2.get_id()
        assert np.all(np.equal(
            edge1.get_position(),
            edge2.get_position(),
        ))
        assert edge1.get_periodicity() == edge2.get_periodicity()
