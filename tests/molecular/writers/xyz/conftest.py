import pytest

import stk

from .case_data import CaseData


@pytest.fixture(
    params=(
        lambda: CaseData(
            molecule=stk.BuildingBlock("BrCCBr", [stk.BromoFactory()]),
            writer=stk.XyzWriter(),
            string=(
                "8\n\nBr -1.379461 1.585348 0.276855\nC -0.729375 -0.231170 "
                "0.138562\nC 0.709973 -0.142141 -0.347528\nBr 1.633541 "
                "0.894493 1.020494\nH -1.386127 -0.819696 -0.517798\nH -0.764"
                "640 -0.635286 1.185943\nH 0.746592 0.491188 -1.260450\nH 1."
                "169497 -1.142737 -0.496078\n"
            ),
        ),
    ),
)
def case_data(request) -> CaseData:
    return request.param()
