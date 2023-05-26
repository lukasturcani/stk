import pytest
import stk

from ..case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda: CaseData(
            factory=stk.RingAmineFactory(),
            molecule=stk.BuildingBlock("NCC(Br)c1c(Br)cccc1"),
            functional_groups=(
                stk.RingAmine(
                    nitrogen=stk.N(0),
                    carbon1=stk.C(1),
                    carbon2=stk.C(2),
                    carbon3=stk.C(4),
                    hydrogen1=stk.H(11),
                    hydrogen2=stk.H(12),
                    hydrogen3=stk.H(15),
                ),
            ),
        ),
    ),
)
def case_data(request) -> CaseData:
    """
    A :class:`.CaseData` instance.

    """

    return request.param()
