import numpy as np
import pytest


@pytest.fixture(
    params=[
        -np.pi / 2,
        np.pi / 2,
    ],
)
def angle(request):
    return request.param


@pytest.fixture(
    params=[
        np.array([0, 1, 0]),
        np.array([1, 0, 0]),
        np.array([1 / np.sqrt(3), 1 / np.sqrt(3), 1 / np.sqrt(3)]),
    ],
)
def axis(request):
    return np.array(request.param)
