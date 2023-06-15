import pytest
import stk

from ....case_data import CaseData
from .utilities import get_closest_point, get_edges, get_fg_position


@pytest.fixture
def tail_3(position, building_block_2):
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

    vertex = stk.polymer.TailVertex(0, position, True)
    return CaseData(
        vertex=vertex,
        edges=(tuple(get_edges(vertex))[0],),
        building_block=building_block_2,
        position=position,
        alignment_tests={
            get_fg0_point: point2,
            get_fg1_point: point1,
        },
        functional_group_edges={1: 0},
    )
