import numpy as np
import pytest
from pytest_lazyfixture import lazy_fixture

# Fixtures must be visible for lazy_fixture() calls.
from .fixtures import *  # noqa


@pytest.fixture(
    params=(
        lazy_fixture("cof1"),
        lazy_fixture("cof2"),
        lazy_fixture("cof3"),
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
