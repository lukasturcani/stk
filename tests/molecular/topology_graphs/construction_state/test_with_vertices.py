from __future__ import annotations

import pytest
import stk

from .utilities import is_clone


@pytest.mark.skip
def test_with_vertices(construction_state, vertices):
    clone = construction_state.clone()
    _test_with_vertices(construction_state, vertices)
    # Test immutability.
    is_clone(construction_state, clone)


def _test_with_vertices(construction_state, vertices):
    clone = construction_state.with_vertices(vertices)
    assert clone.get_num_vertices() == len(vertices)
    for vertex_id in range(len(vertices)):
        assert clone.get_vertex(vertex_id) is vertices[vertex_id]


@pytest.fixture(
    scope="session",
    params=(
        (
            lambda: stk.Vertex(0, [0, 0, 0]),
            lambda: stk.Vertex(1, [10, 0, 0]),
            lambda: stk.Vertex(2, [20, 0, 0]),
        ),
    ),
)
def vertices(request) -> tuple[stk.Vertex, ...]:
    return request.param()
