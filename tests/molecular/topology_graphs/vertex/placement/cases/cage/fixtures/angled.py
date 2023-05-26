from __future__ import annotations

from functools import partial

import numpy as np
import pytest
import stk
from scipy.spatial.distance import euclidean

from ....case_data import CaseData


@pytest.fixture
def angled(
    position: np.ndarray,
    angled_aligner_edge: int,
    building_block_2: stk.BuildingBlock,
) -> CaseData:
    point1, point2 = points = (
        position + [0, 10, 0],
        position + [10, 0, 0],
    )

    def get_fg0_point(building_block: stk.BuildingBlock) -> float:
        return get_closest_point(
            points=points,
            point=get_fg_position(0, building_block),
        )

    def get_fg1_point(building_block: stk.BuildingBlock) -> float:
        return get_closest_point(
            points=points,
            point=get_fg_position(1, building_block),
        )

    vertex = stk.cage.AngledVertex(
        id=0,
        position=position,
        aligner_edge=angled_aligner_edge,
    )

    return CaseData(
        vertex=vertex,
        edges=tuple(get_angled_edges(vertex)),
        building_block=building_block_2,
        position=position,
        alignment_tests={
            get_fg0_point: (point1 if angled_aligner_edge == 0 else point2),
            get_fg1_point: (point2 if angled_aligner_edge == 0 else point1),
        },
        functional_group_edges=(
            {0: 0, 1: 1} if angled_aligner_edge == 0 else {0: 1, 1: 0}
        ),
    )


def get_closest_point(
    points: tuple[np.ndarray, ...],
    point: np.ndarray,
) -> float:
    return min(points, key=partial(euclidean, point))


def get_fg_position(
    id: id,
    building_block: stk.BuildingBlock,
) -> np.ndarray:
    functional_group = next(building_block.get_functional_groups(id))
    return building_block.get_centroid(
        atom_ids=functional_group.get_placer_ids(),
    )


def get_angled_edges(vertex: stk.Vertex) -> stk.Edge:
    vertex2 = stk.Vertex(1, vertex.get_position() + [0, 10, 0])
    vertex3 = stk.Vertex(2, vertex.get_position() + [10, 0, 0])
    yield stk.Edge(0, vertex, vertex2)
    yield stk.Edge(1, vertex, vertex3)


@pytest.fixture(params=(0, 1))
def angled_aligner_edge(request) -> int:
    return request.param


@pytest.fixture(
    params=(
        lambda: stk.BuildingBlock(
            smiles="BrCCNCBr",
            functional_groups=[stk.BromoFactory()],
        ),
    ),
)
def building_block_2(request) -> stk.BuildingBlock:
    return request.param()


@pytest.fixture(
    params=([1, 2, -20],),
)
def position(request) -> np.ndarray:
    """
    The `position` of a vertex.

    """

    return np.array(request.param, dtype=np.float64)
