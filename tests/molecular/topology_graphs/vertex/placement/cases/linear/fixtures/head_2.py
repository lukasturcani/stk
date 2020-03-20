import pytest
import stk

from .utilities import get_closest_point, get_fg_position, get_edges
from ....case_data import CaseData

vertices = stk.molecular.topology_graphs.polymer.linear


@pytest.fixture
def head_2(position, building_block_2):
    point1, point2 = points = (
        position + [-10, 0, 0],
        position + [10, 0, 0],
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

    vertex = vertices._HeadVertex(0, position, False)
    return CaseData(
        vertex=vertex,
        edges=(tuple(get_edges(vertex))[1], ),
        building_block=building_block_2,
        position=position,
        alignment_tests={
            get_fg0_point: point1,
            get_fg1_point: point2,
        },
        functional_group_edges={1: 1},
    )
