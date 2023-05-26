from functools import partial

import pytest
import stk
from scipy.spatial.distance import euclidean

from ....case_data import CaseData


@pytest.fixture
def monodentate(position, building_block_1):
    point1, point2 = points = (
        position + [-10, 0, 0],
        position + [10, 0, 0],
    )

    def get_fg0_point(building_block):
        return get_closest_point(
            points=points,
            point=get_fg_position(0, building_block),
        )

    def get_core_point(building_block):
        return get_closest_point(
            points=points,
            point=get_core_position(building_block),
        )

    vertex = stk.metal_complex.MonoDentateLigandVertex(
        id=0,
        position=position,
    )

    return CaseData(
        vertex=vertex,
        edges=tuple(get_edge(vertex)),
        building_block=building_block_1,
        position=position,
        alignment_tests={
            get_fg0_point: point1,
            get_core_point: point2,
        },
        functional_group_edges=({0: 0}),
    )


def get_closest_point(points, point):
    return min(points, key=partial(euclidean, point))


def get_fg_position(id, building_block):
    functional_group = next(building_block.get_functional_groups(id))
    return building_block.get_centroid(
        atom_ids=functional_group.get_placer_ids(),
    )


def get_core_position(building_block):
    return building_block.get_centroid(
        atom_ids=building_block.get_core_atom_ids(),
    )


def get_edge(vertex):
    vertex2 = stk.Vertex(1, vertex.get_position() + [-1, 0, 0])
    yield stk.Edge(0, vertex, vertex2)


@pytest.fixture(
    params=(
        lambda: stk.BuildingBlock(
            smiles="c1cc2c(cn1)CCCCC2",
            functional_groups=[
                stk.SmartsFunctionalGroupFactory(
                    smarts="[#6]~[#7X2]~[#6]",
                    bonders=(1,),
                    deleters=(),
                ),
            ],
        ),
    ),
)
def building_block_1(request) -> stk.BuildingBlock:
    return request.param()
