import numpy as np
import pytest
import stk

from ...case_data import CaseData


@pytest.fixture(
    params=(
        lambda: CaseData(
            vertex=stk.macrocycle.CycleVertex(
                id=0,
                position=(1, 2, 3),
                flip=True,
                angle=np.pi,
            ),
            id=0,
            position=np.array([1, 2, 3], dtype=np.float64),
            cell=np.array([0, 0, 0]),
        ),
    ),
)
def case_data(request) -> CaseData:
    return request.param()
