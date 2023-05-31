from functools import partial

import pytest
import stk
from scipy.spatial.distance import euclidean

from ....case_data import CaseData


@pytest.fixture
def bidentate(position, building_block_2):
    point1, point2 = points = (
        position + [-10, 0, 0],
        position + [0, -10, 0],
    )

    def get_fg0_point(building_block):
        return get_closest_point(
            points=points,
            point=get_fg_position(0, building_block),
        )

    def get_fg1_point(building_block):
        return get_closest_point(
            points=points,
            point=get_fg_position(1, building_block),
        )

    vertex = stk.metal_complex.BiDentateLigandVertex(
        id=0,
        position=position,
    )

    return CaseData(
        vertex=vertex,
        edges=tuple(get_edges(vertex)),
        building_block=building_block_2,
        position=position,
        alignment_tests={
            get_fg0_point: point1,
            get_fg1_point: point2,
        },
        functional_group_edges=({0: 0, 1: 1}),
    )


def get_closest_point(points, point):
    return min(points, key=partial(euclidean, point))


def get_fg_position(id, building_block):
    functional_group = next(building_block.get_functional_groups(id))
    return building_block.get_centroid(
        atom_ids=functional_group.get_placer_ids(),
    )


def get_edges(vertex):
    vertex2 = stk.Vertex(1, vertex.get_position() + [-1, -1, 0])
    yield stk.Edge(
        id=0,
        vertex1=vertex,
        vertex2=vertex2,
        position=vertex.get_position() + [-1, 0, 0],
    )
    yield stk.Edge(
        id=1,
        vertex1=vertex,
        vertex2=vertex2,
        position=vertex.get_position() + [0, -1, 0],
    )


@pytest.fixture(
    params=(
        lambda: stk.BuildingBlock(
            smiles="C=NC/C=N/Br",
            functional_groups=[
                stk.SmartsFunctionalGroupFactory(
                    smarts="[#6]~[#7X2]~[#35]",
                    bonders=(1,),
                    deleters=(),
                ),
                stk.SmartsFunctionalGroupFactory(
                    smarts="[#6]~[#7X2]~[#6]",
                    bonders=(1,),
                    deleters=(),
                ),
            ],
        ),
    ),
)
def building_block_2(request) -> stk.BuildingBlock:
    return request.param()
