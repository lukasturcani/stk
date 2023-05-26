import pytest
import stk

from .case_data import CaseData


@pytest.fixture(
    params=(
        lambda: CaseData(
            molecule=stk.BuildingBlock("NCCN"),
            bonds=(
                stk.Bond(stk.N(0), stk.C(1), 1),
                stk.Bond(stk.C(1), stk.C(2), 1),
                stk.Bond(stk.C(2), stk.N(3), 1),
                stk.Bond(stk.N(0), stk.H(4), 1),
                stk.Bond(stk.N(0), stk.H(5), 1),
                stk.Bond(stk.C(1), stk.H(6), 1),
                stk.Bond(stk.C(1), stk.H(7), 1),
                stk.Bond(stk.C(2), stk.H(8), 1),
                stk.Bond(stk.C(2), stk.H(9), 1),
                stk.Bond(stk.N(3), stk.H(10), 1),
                stk.Bond(stk.N(3), stk.H(11), 1),
            ),
        ),
    ),
)
def case_data(request) -> CaseData:
    return request.param()
