import pytest
import stk

from .case_data import CaseData


@pytest.fixture(
    params=(
        lambda: CaseData(
            molecule=stk.BuildingBlock("BrCCBr", [stk.BromoFactory()]),
            writer=stk.XyzWriter(),
            string=(
                "8\n\nBr -1.423838 1.561473 0.322335\nC -0.740543 -0.2"
                "57311 0.127980\nC 0.714791 -0.115704 -0.338259\nBr 1."
                "626726 0.889555 1.068701\nH -1.351758 -0.807456 -0.59"
                "3854\nH -0.776931 -0.696380 1.144036\nH 0.769475 0.52"
                "7986 -1.238698\nH 1.182078 -1.102163 -0.492240\n"
            ),
        ),
    ),
)
def case_data(request) -> CaseData:
    return request.param()
