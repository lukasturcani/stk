import pytest
import stk

from .case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            molecule=stk.BuildingBlock('NCCN'),
            path='NCCN.mae',
        ),

        CaseData(
            molecule=stk.BuildingBlock('NCCN'),
            path='NCCN.xyz',
        ),

        CaseData(
            molecule=stk.BuildingBlock('NCCN'),
            path='NCCN.pdb',
        ),

        CaseData(
            molecule=stk.BuildingBlock('NCCN'),
            path='NCCN_with_cell.pdb',
        ),

        CaseData(
            molecule=stk.BuildingBlock('NCCN'),
            path='NCCN.mol',
        ),

        CaseData(
            molecule=stk.BuildingBlock('NCCN'),
            path='NCCN.coord',
        ),

        CaseData(
            molecule=stk.BuildingBlock('NCCN'),
            path='NCCN_with_cell1.coord',
        ),

        CaseData(
            molecule=stk.BuildingBlock('NCCN'),
            path='NCCN_with_cell2.coord',
        ),
    ),
)
def case_data(request):
    """
    A :class:`.CaseData` instance.

    """

    return request.param
