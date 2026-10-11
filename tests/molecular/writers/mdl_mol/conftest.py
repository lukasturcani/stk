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
                " V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 8 7 0 0 0\nM  V30 "
                "BEGIN ATOM\nM  V30 1 Br -1.3795 1.5853 0.2769 0\nM  V30 2 C "
                "-0.7294 -0.2312 0.1386 0\nM  V30 3 C 0.7100 -0.1421 -0.3475 "
                "0\nM  V30 4 Br 1.6335 0.8945 1.0205 0\nM  V30 5 H -1.3861 -0."
                "8197 -0.5178 0\nM  V30 6 H -0.7646 -0.6353 1.1859 0\nM  V30 7"
                " H 0.7466 0.4912 -1.2604 0\nM  V30 8 H 1.1695 -1.1427 -0.4961"
                " 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2\nM  V3"
                "0 2 1 2 3\nM  V30 3 1 3 4\nM  V30 4 1 2 5\nM  V30 5 1 2 6\nM "
                " V30 6 1 3 7\nM  V30 7 1 3 8\nM  V30 END BOND\nM  V30 END CTA"
                "B\nM  END\n\n$$$$\n"
            ),
        ),
    ),
)
def case_data(request) -> CaseData:
    return request.param()
