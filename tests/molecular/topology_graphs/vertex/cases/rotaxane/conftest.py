import numpy as np
import pytest
import stk

from ...case_data import CaseData


@pytest.fixture(
    params=(
        lambda: CaseData(
            vertex=stk.rotaxane.AxleVertex(0, (1, 2, 3)),
            id=0,
            position=np.array([1, 2, 3], dtype=np.float64),
            cell=np.array([0, 0, 0]),
        ),
        lambda: CaseData(
            vertex=stk.rotaxane.CycleVertex(0, (3, 4, 5), True),
            id=0,
            position=np.array([3, 4, 5], dtype=np.float64),
            cell=np.array([0, 0, 0]),
        ),
    ),
)
def case_data(request) -> CaseData:
    return request.param()
