from functools import partial

import numpy as np
import pytest
import stk
from pytest_lazyfixture import lazy_fixture
from scipy.spatial.distance import euclidean

from ....case_data import CaseData


@pytest.fixture(
    params=(
        lazy_fixture("nonlinear_3"),
        lazy_fixture("nonlinear_4"),
    ),
)
def nonlinear(request):
    return request.param


@pytest.fixture
def nonlinear_3(position, aligner_3, building_block_3):
    return _nonlinear(position, aligner_3, building_block_3)


@pytest.fixture(params=(0, 1, 2))
def aligner_3(request):
    """
    The `aligner_edge` parameter for vertices with 3 edges.

    """

    return request.param


@pytest.fixture(
    params=(
        lambda: stk.BuildingBlock(
            smiles="BrCC(Br)C(CCCC)NCBr",
            functional_groups=[stk.BromoFactory()],
        ),
    ),
)
def building_block_3(request) -> stk.BuildingBlock:
    """
    A :class:`.BuildingBlock` with 3 functional groups.

    """

    return request.param()


@pytest.fixture
def nonlinear_4(position, aligner_4, building_block_4):
    return _nonlinear(position, aligner_4, building_block_4)


@pytest.fixture(params=(0, 1, 2, 3))
def aligner_4(request):
    """
    The `aligner_edge` parameter for vertices with 4 edges.

    """

    return request.param


@pytest.fixture(
    params=(
        lambda: stk.BuildingBlock(
            smiles="BrC(CCCCCC)C(Br)CC(Br)CNCBr",
            functional_groups=[stk.BromoFactory()],
        ),
    ),
)
def building_block_4(request) -> stk.BuildingBlock:
    """
    A :class:`.BuildingBlock` with 4 functional groups.

    """

    return request.param()


@pytest.fixture(
    params=(
        [1, 2, -20],
        [1, 2, 20],
    ),
)
def position(request):
    """
    The `position` of a vertex.

    """

    return np.array(request.param, dtype=np.float64)


def _nonlinear(position, aligner_edge, building_block):
    """
    Return a test case for a nonlinear COF vertex.

    Parameters
    ----------
    position : :class:`numpy.ndarray`
        The position of the vertex.

    aligner_edge : :class:`int`
        The aligner edge of the vertex.

    building_block : :class:`.BuildingBlock`
        The building block placed on the vertex.

    Returns
    -------
    :class:`.CaseData`
        The test case.

    """

    building_block = order_functional_groups(building_block)
    n = building_block.get_num_functional_groups()
    points = tuple(get_points(position, n))
    alignment_tests = {
        partial(get_fg_point, points, 0): points[aligner_edge],
        get_normal: (
            np.array([0, 0, 1]) if position[2] > 0 else np.array([0, 0, -1])
        ),
    }

    vertex = stk.cage.NonLinearVertex(
        id=0,
        position=position,
        aligner_edge=aligner_edge,
    )

    return CaseData(
        vertex=vertex,
        edges=tuple(get_nonlinear_edges(n, vertex)),
        building_block=building_block,
        position=position,
        alignment_tests=alignment_tests,
        functional_group_edges={
            fg_id: (fg_id + aligner_edge) % n for fg_id in range(n)
        },
    )


def get_points(center, num_points):
    """
    Yield equally spaced points on a circle of radius 10.

    Parameters
    ----------
    center : :class:`numpy.ndarray`
        The center of circle.

    num_points : :class:`int`
        The number of points to make.

    Yields
    ------
    :class:`numpy.ndarray`
        A point on a circle.

    """

    # Generate points in counter clockwise order if the building block
    # is to be placed upside down. This keeps the relative ordering
    # of points the same from the perspective of the building block.
    direction = 1 if center[2] > 0 else -1

    # Take slice to account for case where rounding errors cause
    # extra theta.
    thetas = np.arange(
        start=0,
        stop=direction * 2 * np.pi,
        step=direction * 2 * np.pi / num_points,
    )[:num_points]
    for theta in thetas:
        yield center + 10 * np.array([np.sin(theta), np.cos(theta), 0.0])


def get_fg_point(points, fg_id, building_block):
    """
    Get the point in `points` closest to a functional group.

    Parameters
    ----------
    points : :class:`iterable` of :class:`numpy.ndarray`
        The points from which the closest one is picked.

    fg_id : :class:`int`
        The id of a functional group of `building_block` whose
        distance to each point in `points` is evaluated.

    building_block : :class:`.BuildingBlock`
        The building block which owns the functional group.

    Returns
    -------
    :class:`numpy.ndarray`
        The point in `points` closest to functional group `fg_id`.

    """

    fg = next(building_block.get_functional_groups(fg_id))
    fg_position = building_block.get_centroid(
        atom_ids=fg.get_placer_ids(),
    )
    return min(points, key=partial(euclidean, fg_position))


def get_normal(building_block):
    """
    Get a normal to a plane crossing the building block.

    The normal is defined by the plane of best fit to the placer atoms
    of the building block. It is defined such that the it always forms
    an acute angle with the vector running from the centroid of the
    placer atoms to the centroid of the building block.

    Parameters
    ----------
    building_block : :.BuildingBlock`
        The building block whose normal should be calculated.

    Returns
    -------
    :class:`numpy.ndarray`
        The normal of `building_block`.

    """

    placer_centroid = building_block.get_centroid(
        atom_ids=building_block.get_placer_ids(),
    )
    normal = stk.get_acute_vector(
        reference=building_block.get_centroid() - placer_centroid,
        vector=building_block.get_plane_normal(
            atom_ids=building_block.get_placer_ids(),
        ),
    )
    if np.allclose(normal, [0, 0, 1], atol=1e-13):
        return np.array([0, 0, 1])
    if np.allclose(normal, [0, 0, -1], atol=1e-13):
        return np.array([0, 0, -1])
    return normal


def get_nonlinear_edges(num_edges, vertex):
    """
    Yield edges placed in a circle around `vertex`.

    Parameters
    ----------
    num_edges : :class:`int`
        The number of edges to yield.

    vertex : :class:`.Vertex`
        The vertex which needs edges.

    Yields
    ------
    :class:`.Edge`
        An edge connected to `vertex`.

    """

    for id_, point in enumerate(
        get_points(center=vertex.get_position(), num_points=num_edges)
    ):
        yield stk.Edge(id_, vertex, stk.Vertex(id_ + 1, point))


def order_functional_groups(building_block):
    """
    Get a building block with ordered functional groups.

    Parameters
    ----------
    building_block : :class:`.BuildingBlock`
        The building block whose functional groups should be
        ordered.

    Returns
    -------
    :class:`.BuildingBlock`
        A building block with ordered functional groups.

    """

    return building_block.with_functional_groups(
        functional_groups=sorted(
            building_block.get_functional_groups(),
            key=functional_group_angle(building_block),
        ),
    )


def functional_group_angle(building_block):
    fg0 = next(building_block.get_functional_groups(0))
    fg0_position = building_block.get_centroid(
        atom_ids=fg0.get_placer_ids(),
    )
    placer_centroid = building_block.get_centroid(
        atom_ids=building_block.get_placer_ids(),
    )
    fg0_direction = fg0_position - placer_centroid
    centroid = building_block.get_centroid(
        atom_ids=building_block.get_placer_ids(),
    )
    normal = stk.get_acute_vector(
        reference=building_block.get_centroid() - placer_centroid,
        vector=building_block.get_plane_normal(
            atom_ids=building_block.get_placer_ids(),
        ),
    )
    axis = -np.cross(normal, fg0_direction)

    def inner(functional_group):
        position = building_block.get_centroid(
            atom_ids=functional_group.get_placer_ids(),
        )
        fg_direction = position - centroid
        theta = stk.vector_angle(fg0_direction, fg_direction)

        projection = fg_direction @ axis
        if theta > 0 and projection < 0:
            return 2 * np.pi - theta
        return theta

    return inner
