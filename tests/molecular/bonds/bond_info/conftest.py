import pytest
import stk

from .case_data import CaseData


@pytest.fixture(
    params=(stk.Bond(stk.C(0), stk.N(1), 1),),
)
def bond(request):
    """
    A :class:`.Bond` instance.

    """

    return request.param


@pytest.fixture(
    params=(
        stk.BuildingBlock("NCCN"),
        None,
    ),
)
def building_block(request):
    """
    A valid building block for a :class:`.BondInfo`.

    """

    return request.param


@pytest.fixture(
    params=(
        1,
        None,
    ),
)
def building_block_id(request):
    """
    A valid building block id.

    """

    return request.param


@pytest.fixture
def case_data(bond, building_block, building_block_id):
    """
    A :class:`.CaseData` instance.

    """

    return CaseData(bond, building_block, building_block_id)
