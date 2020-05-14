import pytest
import numpy as np
import stk

from ...case_data import CaseData

vertices = stk.metal_complex.vertices


@pytest.fixture(
    params=(
        CaseData(
            vertex=vertices._MetalVertex(0, (1, 2, 3)),
            id=0,
            position=np.array([1, 2, 3], dtype=np.float64),
            cell=np.array([0, 0, 0]),
        ),
        CaseData(
            vertex=vertices._MonoDentateLigandVertex(0, (1, 2, 3)),
            id=0,
            position=np.array([1, 2, 3], dtype=np.float64),
            cell=np.array([0, 0, 0]),
        ),
        CaseData(
            vertex=vertices._BiDentateLigandVertex(0, (1, 2, 3)),
            id=0,
            position=np.array([1, 2, 3], dtype=np.float64),
            cell=np.array([0, 0, 0]),
        ),
    ),
)
def case_data(request):
    return request.param
