import pytest
import stk

from .case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            molecule=stk.BuildingBlock('NCCN'),
            atoms=(
                stk.N(0),
                stk.C(1),
                stk.C(2),
                stk.N(3),
                stk.H(4),
                stk.H(5),
                stk.H(6),
                stk.H(7),
                stk.H(8),
                stk.H(9),
                stk.H(10),
                stk.H(11),
            ),
        ),
    ),
)
def case_data(request):
    """
    A :class:`.CaseData` instance.

    """

    return request.param
