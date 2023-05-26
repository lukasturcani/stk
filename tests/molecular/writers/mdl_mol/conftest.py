import pytest
import stk

from .case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda: CaseData(
            molecule=stk.BuildingBlock("BrCCBr", [stk.BromoFactory()]),
            writer=stk.MolWriter(),
            string=(
                "\n     RDKit          3D\n\n  0  0  0  0  0  0  0  0 "
                " 0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 8 7 "
                "0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 Br -1.4238 1.5615 "
                "0.3223 0\nM  V30 2 C -0.7405 -0.2573 0.1280 0\nM  V30"
                " 3 C 0.7148 -0.1157 -0.3383 0\nM  V30 4 Br 1.6267 0.8"
                "896 1.0687 0\nM  V30 5 H -1.3518 -0.8075 -0.5939 0\nM"
                "  V30 6 H -0.7769 -0.6964 1.1440 0\nM  V30 7 H 0.7695"
                " 0.5280 -1.2387 0\nM  V30 8 H 1.1821 -1.1022 -0.4922 "
                "0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2"
                "\nM  V30 2 1 2 3\nM  V30 3 1 3 4\nM  V30 4 1 2 5\nM  "
                "V30 5 1 2 6\nM  V30 6 1 3 7\nM  V30 7 1 3 8\nM  V30 E"
                "ND BOND\nM  V30 END CTAB\nM  END\n\n$$$$\n"
            ),
        ),
    ),
)
def case_data(request) -> CaseData:
    return request.param()
