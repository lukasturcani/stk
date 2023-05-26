import numpy as np
import pytest
import stk
from pytest_lazyfixture import lazy_fixture

# Fixtures need to be visible for lazy_fixture() calls.
from .fixtures import *  # noqa


@pytest.fixture(
    params=(
        lazy_fixture("flip"),
        lazy_fixture("no_flip"),
    ),
)
def case_data(request):
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
        np.pi / 2,
        0,
        np.pi,
    ),
)
def angle(request):
    return request.param
