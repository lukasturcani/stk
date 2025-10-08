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
                "\n     RDKit          3D\n\n  0  0  0  0  0  0  0  0  0  0999"
                " V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 8 7 0 0 0\nM  V30 BE"
                "GIN ATOM\nM  V30 1 Br -1.4302 1.5499 0.3629 0\nM  V30 2 C -0."
                "7383 -0.2592 0.1157 0\nM  V30 3 C 0.7182 -0.1112 -0.3446 0\nM"
                "  V30 4 Br 1.6324 0.8520 1.0901 0\nM  V30 5 H -1.3636 -0.7969"
                " -0.5986 0\nM  V30 6 H -0.7775 -0.7041 1.1291 0\nM  V30 7 H 0"
                ".7784 0.5649 -1.2208 0\nM  V30 8 H 1.1805 -1.0955 -0.5339 0\n"
                "M  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2\nM  V30 2 "
                "1 2 3\nM  V30 3 1 3 4\nM  V30 4 1 2 5\nM  V30 5 1 2 6\nM  V30"
                " 6 1 3 7\nM  V30 7 1 3 8\nM  V30 END BOND\nM  V30 END CTAB\nM"
                "  END\n\n$$$$\n"
            ),
        ),
    ),
)
def case_data(request) -> CaseData:
    return request.param()
