import pytest
import stk
import numpy as np
from pytest_lazyfixture import lazy_fixture

# Fixtures must by visible for lazy_fixture() calls.
from .fixtures import *  # noqa

vertices = stk.molecular.topology_graphs.polymer.linear


@pytest.fixture(
    params=(
        lazy_fixture('center'),
        lazy_fixture('head_1'),
        lazy_fixture('head_2'),
        lazy_fixture('head_3'),
        lazy_fixture('tail_1'),
        lazy_fixture('tail_2'),
        lazy_fixture('tail_3'),
    ),
)
def case_data(request):
    """
    A :class:`.CaseData` instance.

    """

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
        stk.BuildingBlock(
            smiles='BrCC',
            functional_groups=[stk.BromoFactory()],
        ),
    ),
)
def building_block_1(request):
    """
    A :class:`.BuildingBlock` with 1 functional group.

    """

    return request.param


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
