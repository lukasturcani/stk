import pytest
import numpy as np
import stk
from functools import partial
from scipy.spatial.distance import euclidean


from ..._test_case import _TestCase

vertices = stk.molecular.topology_graphs.cof.vertices


@pytest.fixture
def nonlinear(position, nonlinear_aligner_edge, building_block_3):

    point1, point2, point3 = points = (
        position + [0, 10, 0],
        position + [7, -7, 0],
        position + [-7, -7, 0],
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

    def get_fg2_point(building_block):
        return get_closest_point(
            points=points,
            point=get_fg_position(2, building_block),
        )

    def get_normal(building_block):
        placer_centroid = building_block.get_centroid(
            atom_ids=building_block.get_placer_ids(),
        )
        normal = building_block.get_plane_normal(
            atom_ids=building_block.get_placer_ids(),
        )
        normal = stk.get_acute_vector(
            reference=building_block.get_centroid() - placer_centroid,
            vector=normal,
        )
        if np.allclose(normal, [0, 0, 1], atol=1e-13):
            return np.array([0, 0, 1])
        return normal

    vertex = vertices._NonLinearCofVertex(
        id=0,
        position=position,
        aligner_edge=nonlinear_aligner_edge,
        cell=[0, 0, 0],
    )

    return _TestCase(
        vertex=vertex,
        edges=tuple(get_nonlinear_edges(vertex)),
        building_block=building_block_3,
        position=position,
        points={
            get_fg0_point: get_expected_point(
                fg=0,
                aligner_edge=nonlinear_aligner_edge,
                points=points,
            ),
            get_fg1_point: get_expected_point(
                fg=1,
                aligner_edge=nonlinear_aligner_edge,
                points=points,
            ),
            get_fg2_point: get_expected_point(
                fg=2,
                aligner_edge=nonlinear_aligner_edge,
                points=points,
            ),
            get_normal: np.array([0, 0, 1]),
        },
        functional_group_edges=get_functional_group_edges(
            aligner_edge=nonlinear_aligner_edge,
        ),
    )


def get_functional_group_edges(aligner_edge):
    return {
        0: {0: 0, 1: 1, 2: 2},
        1: {0: 1, 1: 2, 2: 0},
        2: {0: 2, 1: 0, 2: 1},
    }[aligner_edge]


def get_expected_point(fg, aligner_edge, points):
    point1, point2, point3 = points
    return {
        (0, 0): point1,
        (0, 1): point2,
        (0, 2): point3,
        (1, 0): point2,
        (1, 1): point3,
        (1, 2): point1,
        (2, 0): point3,
        (2, 1): point1,
        (2, 2): point2,
    }[(fg, aligner_edge)]


def get_closest_point(points, point):
    return min(points, key=partial(euclidean, point))


def get_fg_position(id, building_block):
    functional_group = next(building_block.get_functional_groups(id))
    return building_block.get_centroid(
        atom_ids=functional_group.get_placer_ids(),
    )


def get_nonlinear_edges(vertex):
    vertex2 = stk.Vertex(1, vertex.get_position() + [0, 10, 0])
    vertex3 = stk.Vertex(2, vertex.get_position() + [7, -7, 0])
    vertex4 = stk.Vertex(3, vertex.get_position() + [7, 7, 0])
    yield stk.Edge(0, vertex, vertex2)
    yield stk.Edge(1, vertex, vertex3)
    yield stk.Edge(2, vertex, vertex4)


@pytest.fixture(params=(0, 1, 2))
def nonlinear_aligner_edge(request):
    return request.param


@pytest.fixture(
    params=(
        stk.BuildingBlock(
            smiles='BrCC(Br)CNCBr',
            functional_groups=[stk.BromoFactory()],
        ),
    ),
)
def building_block_3(request):
    return order_functional_groups(request.param)


def order_functional_groups(building_block):
    building_block = building_block.with_centroid(
        position=[0, 0, 0],
        atom_ids=building_block.get_placer_ids(),
    )
    normal = building_block.get_plane_normal(
        atom_ids=building_block.get_placer_ids(),
    )
    normal = stk.get_acute_vector(
        reference=building_block.get_centroid(),
        vector=normal,
    )
    building_block = building_block.with_rotation_between_vectors(
        start=normal,
        target=[0, 0, 1],
        origin=np.array([0, 0, 0]),
    )
    fg = next(building_block.get_functional_groups(0))
    fg_centroid = building_block.get_centroid(fg.get_placer_ids())
    building_block = (
        building_block.with_rotation_to_minimize_angle(
            start=fg_centroid,
            target=[0, 1, 0],
            axis=[0, 0, 1],
            origin=np.array([0, 0, 0]),
        )
    )
    return building_block.with_functional_groups((
        next(building_block.get_functional_groups(0)),
        max(
            building_block.get_functional_groups(),
            key=get_x_coord(building_block),
        ),
        min(
            building_block.get_functional_groups(),
            key=get_x_coord(building_block),
        ),
    ))


def get_x_coord(building_block):

    def inner(functional_group):
        return building_block.get_centroid(
            atom_ids=functional_group.get_placer_ids()
        )[0]

    return inner
