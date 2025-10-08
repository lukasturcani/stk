import pytest

import stk

from .case_data import CaseData


@pytest.fixture(
    params=(
        lambda: CaseData(
            molecule=stk.BuildingBlock("BrCCBr", [stk.BromoFactory()]),
            writer=stk.XyzWriter(),
            string=(
                "8\n\nBr -1.430189 1.549908 0.362949\nC -0.738254 -0.259169 0."
                "115667\nC 0.718233 -0.111206 -0.344580\nBr 1.632362 0.852014 "
                "1.090150\nH -1.363556 -0.796869 -0.598581\nH -0.777513 -0.704"
                "050 1.129083\nH 0.778398 0.564873 -1.220778\nH 1.180518 -1.09"
                "5500 -0.533910\n"
            ),
        ),
    ),
)
def case_data(request) -> CaseData:
    return request.param()
