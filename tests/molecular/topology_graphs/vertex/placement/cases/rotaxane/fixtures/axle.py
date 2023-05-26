import numpy as np
import pytest
import stk

from ....case_data import CaseData


@pytest.fixture(
    params=(
        lambda: CaseData(
            vertex=stk.rotaxane.AxleVertex(0, (1, 2, 3)),
            edges=(),
            building_block=stk.BuildingBlock("BrCCBr"),
            position=np.array([1, 2, 3], dtype=np.float64),
            alignment_tests={},
            functional_group_edges={},
        ),
    ),
)
def axle(request) -> CaseData:
    return request.param()
