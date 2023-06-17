import pytest
import stk

from ....case_data import CaseData
from .utilities import get_centroid, get_closest_point, get_edges


@pytest.fixture
def tail_1(position, flip, building_block_1):
    point1, point2 = points = (
        position + [-10, 0, 0],
        position + [10, 0, 0],
    )

    def get_centroid_point(building_block):
        return get_closest_point(
            points=points,
            point=get_centroid(building_block),
        )

    vertex = stk.polymer.TailVertex(0, position, flip)
    return CaseData(
        vertex=vertex,
        edges=(tuple(get_edges(vertex))[0],),
        building_block=building_block_1,
        position=position,
        alignment_tests={get_centroid_point: point2},
        functional_group_edges={0: 0},
    )
