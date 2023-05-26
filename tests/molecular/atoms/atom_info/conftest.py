import pytest
import stk

from .case_data import CaseData


@pytest.fixture(
    params=(stk.C(0),),
)
def atom(request):
    """
    A :class:`.Atom` instance.

    """

    return request.param


@pytest.fixture(
    params=(
        stk.C(0),
        None,
    ),
)
def building_block_atom(request):
    """
    A valid building block atom for :class:`.AtomInfo`.

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
    A valid building block for :class:`.AtomInfo`.

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
    A building block id.

    """

    return request.param


@pytest.fixture
def case_data(
    atom,
    building_block_atom,
    building_block,
    building_block_id,
):
    """
    A :class:`.CaseData` instance.

    """

    return CaseData(
        atom=atom,
        building_block_atom=building_block_atom,
        building_block=building_block,
        building_block_id=building_block_id,
    )
