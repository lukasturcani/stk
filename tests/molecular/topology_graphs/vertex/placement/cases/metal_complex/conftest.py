import numpy as np
import pytest
from pytest_lazyfixture import lazy_fixture

# Fixtures need to be visible for lazy_fixture() calls.
from .fixtures import *  # noqa


@pytest.fixture(
    params=(
        lazy_fixture("metal"),
        lazy_fixture("monodentate"),
        lazy_fixture("bidentate"),
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
