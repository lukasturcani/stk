import pytest
import numpy as np
import stk

from ....case_data import CaseData

vertices = stk.molecular.topology_graphs.host_guest.vertices


@pytest.fixture(
    params=(
        CaseData(
            vertex=vertices._HostVertex(
                id=0,
                position=(1, 2, 3),
            ),
            edges=(),
            building_block=stk.BuildingBlock('NCCN'),
            position=np.array([1, 2, 3], dtype=np.float64),
            alignment_tests={},
            functional_groups_edges={},
        ),
    ),
)
def host(request):
    return request.param
