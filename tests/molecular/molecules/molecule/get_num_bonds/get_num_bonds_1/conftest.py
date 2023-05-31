import pytest
import stk

from .case_data import CaseData


@pytest.fixture(
    params=(
        lambda: CaseData(molecule=stk.BuildingBlock("NCCN"), num_bonds=11),
    ),
)
def case_data(request) -> CaseData:
    """
    A :class:`.CaseData` instance.

    """

    return request.param()
