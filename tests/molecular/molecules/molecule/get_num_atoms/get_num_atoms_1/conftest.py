import pytest
import stk

from .case_data import CaseData


@pytest.fixture(
    params=(
        lambda: CaseData(
            molecule=stk.BuildingBlock("NCCN"),
            num_atoms=12,
        ),
    ),
)
def case_data(request) -> CaseData:
    return request.param()
