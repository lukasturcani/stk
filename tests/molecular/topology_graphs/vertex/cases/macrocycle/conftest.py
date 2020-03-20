import pytest
import numpy as np
from stk.macrocycle import vertices

from ...case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            vertex=vertices._CycleVertex(0, (1, 2, 3), True, np.pi),
            id=0,
            position=np.array([1, 2, 3], dtype=np.float64),
            cell=np.array([0, 0, 0]),
        ),
    ),
)
def case_data(request):
    return request.param
