import numpy as np
import pytest
from pytest_lazy_fixtures import lf

# Fixtures must be visible for lazy_fixture() calls.
from .fixtures import *


@pytest.fixture(
    params=(
        lf("cof1"),
        lf("cof2"),
        lf("cof3"),
    ),
)
def case_data(request):
    return request.param


@pytest.fixture(params=(0,))
def aligner_edge(request):
    return request.param


@pytest.fixture(params=(0,))
def id(request):
    return request.param


@pytest.fixture(
    params=(
        np.array([0, 0, 0]),
        np.array([-20, 1, 21]),
    ),
)
def cell(request):
    return np.array(request.param)


@pytest.fixture(
    params=(
        [0, 0, 0],
        [1, 2, -20],
    ),
)
def position(request):
    return np.array(request.param, dtype=np.float64)
