import numpy as np
import pytest
from pytest_lazy_fixtures import lf

import stk

# Fixtures must by visible for lazy_fixture() calls.
from .fixtures import *


@pytest.fixture(
    params=(
        lf("center"),
        lf("head_1"),
        lf("head_2"),
        lf("head_3"),
        lf("tail_1"),
        lf("tail_2"),
        lf("tail_3"),
    ),
)
def case_data(request):
    """
    A :class:`.CaseData` instance.

    """

    return request.param


@pytest.fixture(
    params=(
        lambda: stk.BuildingBlock(
            smiles="BrCCNCBr",
            functional_groups=[stk.BromoFactory()],
        ),
    ),
)
def building_block_2(request) -> stk.BuildingBlock:
    """
    A :class:`.BuildingBlock` with 2 functional groups.

    """

    return request.param()


@pytest.fixture(
    params=(
        lambda: stk.BuildingBlock(
            smiles="BrCC",
            functional_groups=[stk.BromoFactory()],
        ),
    ),
)
def building_block_1(request) -> stk.BuildingBlock:
    """
    A :class:`.BuildingBlock` with 1 functional group.

    """

    return request.param()


@pytest.fixture(params=(True, False))
def flip(request):
    return request.param


@pytest.fixture(
    params=(
        [0, 0, 0],
        [1, 2, -20],
    ),
)
def position(request):
    return np.array(request.param, dtype=np.float64)
