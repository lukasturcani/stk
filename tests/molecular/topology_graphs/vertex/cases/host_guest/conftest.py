import numpy as np
import pytest
import stk

from ...case_data import CaseData


@pytest.fixture(
    params=(
        lambda: CaseData(
            vertex=stk.host_guest.HostVertex(0, (1, 2, 3)),
            id=0,
            position=np.array([1, 2, 3], dtype=np.float64),
            cell=np.array([0, 0, 0]),
        ),
        lambda: CaseData(
            vertex=stk.host_guest.GuestVertex(
                id=0,
                position=(1, 2, 3),
                start=(4, 5, 6),
                target=(7, 8, 0),
            ),
            id=0,
            position=np.array([1, 2, 3], dtype=np.float64),
            cell=np.array([0, 0, 0]),
        ),
    ),
)
def case_data(request) -> CaseData:
    return request.param()
