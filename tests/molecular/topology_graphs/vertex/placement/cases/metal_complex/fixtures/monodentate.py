import pytest
import numpy as np
import stk
from functools import partial
from scipy.spatial.distance import euclidean

from ....case_data import CaseData

vertices = stk.metal_complex.vertices


@pytest.fixture
def monodentate(position, building_block_1):

    point1 = points = (position + [-10, 0, 0], )

    def get_fg0_point(building_block):
        return get_closest_point(
            points=points,
            point=get_fg_position(0, building_block),
        )

    vertex = vertices._MonoDentateLigandVertex(
        id=0,
        position=position,
    )

    return CaseData(
        vertex=vertex,
        edges=tuple(get_edge(vertex)),
        building_block=building_block_1,
        position=position,
        alignment_tests={get_fg0_point: point1},
        functional_group_edges=({0: 0}),
    )


def get_closest_point(points, point):
    return min(points, key=partial(euclidean, point))


def get_fg_position(id, building_block):
    functional_group = next(building_block.get_functional_groups(id))
    return building_block.get_centroid(
        atom_ids=functional_group.get_placer_ids(),
    )


def get_edge(vertex):
    vertex2 = stk.Vertex(1, vertex.get_position() + [-1, 0, 0])
    yield stk.Edge(0, vertex, vertex2)


@pytest.fixture(
    params=(
        stk.BuildingBlock(
            smiles='c1cc2c(cn1)CCCCC2',
            functional_groups=[
                stk.SmartsFunctionalGroupFactory(
                    smarts='[#6]~[#7X2]~[#6]',
                    bonders=(1, ),
                    deleters=(),
                ),
            ],
        ),
    ),
)
def building_block_1(request):
    return request.param


@pytest.fixture(
    params=(
        [1, 2, -20],
    ),
)
def position(request):
    """
    The `position` of a vertex.

    """

    return np.array(request.param, dtype=np.float64)
