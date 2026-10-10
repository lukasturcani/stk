import numpy as np
import pytest
from pytest_lazy_fixtures import lf

# Fixtures need to be visible for lazy_fixture() calls.
from .fixtures import *


@pytest.fixture(
    params=(
        lf("metal"),
        lf("monodentate"),
        lf("bidentate"),
    ),
)
def case_data(request):
    """
    A :class:`.CaseData` instance.

    """

    return request.param


@pytest.fixture(
    params=([1, 2, -20],),
)
def position(request):
    """
    The `position` of a vertex.

    """

    return np.array(request.param, dtype=np.float64)
