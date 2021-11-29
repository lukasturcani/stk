import numpy as np
import pytest

from .case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(np.array([0, 0, 0]), 2, np.array([0, 0, 0])),
        CaseData(np.array([0, 0, 0]), 0, np.array([0, 0, 0])),
        CaseData(np.array([1, -1, 0]), 0, np.array([0, 0, 0])),
        CaseData(np.array([1, -1, 2]), 2, np.array([2, -2, 4])),
        CaseData(np.array([2, 4, -5]), -3, np.array([-6, -12, 15])),
        CaseData(
            start=np.array([1, -1, 2]),
            scale=np.array([10, 2, -2]),
            target=np.array([10, -2, -4]),
        ),
    ),
)
def case_data(request):
    """
    A :class:`.CaseData` instance.

    """

    return request.param
