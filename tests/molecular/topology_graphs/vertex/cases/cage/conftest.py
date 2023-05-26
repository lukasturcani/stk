from __future__ import annotations

import numpy as np
import pytest
import stk
from pytest_lazyfixture import lazy_fixture

from ...case_data import CaseData


@pytest.fixture(
    params=(
        lazy_fixture("cage1"),
        lazy_fixture("cage2"),
        lazy_fixture("cage3"),
    ),
)
def case_data(request):
    return request.param


@pytest.fixture
def cage1(cls, position):
    return CaseData(
        vertex=cls(
            id=0,
            position=position,
        ),
        id=0,
        position=position,
        cell=np.array([0, 0, 0]),
    )


@pytest.fixture
def cage2(init_at_center, vertices_):
    return CaseData(
        vertex=init_at_center(
            id=0,
            vertices=vertices_,
        ),
        id=0,
        position=(sum(v.get_position() for v in vertices_) / len(vertices_)),
        cell=np.array([0, 0, 0]),
    )


@pytest.fixture
def cage3(position):
    return CaseData(
        vertex=stk.cage.UnaligningVertex(
            id=0,
            position=position,
        ),
        id=0,
        position=position,
        cell=np.array([0, 0, 0]),
    )


@pytest.fixture(
    params=(
        [0, 0, 0],
        [1, 2, -20],
    ),
)
def position(request):
    return np.array(request.param, dtype=np.float64)


@pytest.fixture(
    params=(
        stk.cage.LinearVertex,
        stk.cage.NonLinearVertex,
        stk.cage.UnaligningVertex,
        stk.cage.AngledVertex,
    ),
)
def cls(request) -> stk.Vertex:
    return request.param


@pytest.fixture(
    params=(
        stk.cage.LinearVertex.init_at_center,
        stk.cage.NonLinearVertex.init_at_center,
        stk.cage.UnaligningVertex.init_at_center,
        stk.cage.AngledVertex.init_at_center,
    ),
)
def init_at_center(request) -> stk.Vertex:
    return request.param


@pytest.fixture(
    params=(
        lambda: (stk.Vertex(0, [1, 2, 3]),),
        lambda: (
            stk.Vertex(0, [1, 2, 3]),
            stk.Vertex(1, [-1, 2, -32]),
        ),
    ),
)
def vertices_(request) -> tuple[stk.Vertex, ...]:
    return request.param()
