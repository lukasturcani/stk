import pytest
import numpy as np
from pytest_lazyfixture import lazy_fixture

# Fixtures must be visible for lazy_fixture() calls.
from .fixtures import *  # noqa


@pytest.fixture(
    params=(
        lazy_fixture('linear_placement'),
    ),
)
def placement_test_case(request):
    return request.param


def test_placement(placement_test_case):
    _test_placement(
        vertex=placement_test_case.vertex,
        edges=placement_test_case.edges,
        building_block=placement_test_case.building_block,
        position=placement_test_case.position,
        functional_group_edges=(
            placement_test_case.functional_group_edges
        ),
    )


def _test_placement(
    vertex,
    edges,
    building_block,
    position,
    functional_group_edges,
):
    position_matrix = vertex.place_building_block(building_block)
    building_block = building_block.with_position_matrix(
        position_matrix=position_matrix,
    )
    assert np.allclose(
        a=building_block.get_centroid(building_block.get_placer_ids()),
        b=position,
        atol=1e-14,
    )
    result = vertex.map_functional_groups_to_edges(
        building_block=building_block,
        edges=edges,
    )
    assert result == functional_group_edges
