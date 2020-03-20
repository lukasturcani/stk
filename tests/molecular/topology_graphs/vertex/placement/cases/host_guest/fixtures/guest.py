import pytest
import numpy as np
import stk

from ....case_data import CaseData

vertices = stk.molecular.topology_graphs.host_guest.vertices


def get_aligned_building_block(building_block, target):
    return building_block.with_rotation_between_vectors(
        start=building_block.get_direction(),
        target=target,
        origin=building_block.get_centroid(
            atom_ids=building_block.get_placer_ids(),
        ),
    )


def get_direction(building_block):
    direction = building_block.get_direction()
    if np.allclose(direction, [1, 0, 0], atol=1e-13):
        return np.array([1, 0, 0], dtype=np.float64)
    return direction


@pytest.fixture(
    params=(
        CaseData(
            vertex=vertices._GuestVertex(
                id=0,
                position=(1, 2, 3),
                start=(1, 2, 3),
                target=(1, 0, 0),
            ),
            edges=(),
            building_block=get_aligned_building_block(
                building_block=stk.BuildingBlock('BrCCBr'),
                target=(1, 2, 3),
            ),
            position=np.array([1, 2, 3], dtype=np.float64),
            alignment_tests={
                get_direction: np.array([1, 0, 0], dtype=np.float64),
            },
            functional_group_edges={},
        ),
    ),
)
def guest(request):
    return request.param
