import pytest
import numpy as np
from pytest_lazyfixture import lazy_fixture
import stk

# Fixtures need to be visible for lazy_fixture() calls.
from .fixtures import *  # noqa


@pytest.fixture(
    params=(
        lazy_fixture('flip'),
        lazy_fixture('no_flip'),
    ),
)
def case_data(request):
    return request.param


@pytest.fixture(
    params=(
        stk.BuildingBlock(
            smiles='BrCCNCBr',
            functional_groups=[stk.BromoFactory()],
        ),
    ),
)
def building_block_2(request):
    """
    A :class:`.BuildingBlock` with 2 functional groups.

    """

    return request.param


@pytest.fixture(
    params=(
        np.pi/2,
        0,
        np.pi,
    ),
)
def angle(request):
    return request.param
