import pytest
import stk

from .case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda: CaseData(
            molecule=stk.BuildingBlock("NCCN"),
            path="NCCN.mae",
        ),
        lambda: CaseData(
            molecule=stk.BuildingBlock("NCCN"),
            path="NCCN.xyz",
        ),
        lambda: CaseData(
            molecule=stk.BuildingBlock("NCCN"),
            path="NCCN.pdb",
        ),
        lambda: CaseData(
            molecule=stk.BuildingBlock("NCCN"),
            path="NCCN_with_cell.pdb",
        ),
        lambda: CaseData(
            molecule=stk.BuildingBlock("NCCN"),
            path="NCCN.mol",
        ),
        lambda: CaseData(
            molecule=stk.BuildingBlock("NCCN"),
            path="NCCN.coord",
        ),
        lambda: CaseData(
            molecule=stk.BuildingBlock("NCCN"),
            path="NCCN_with_cell1.coord",
        ),
        lambda: CaseData(
            molecule=stk.BuildingBlock("NCCN"),
            path="NCCN_with_cell2.coord",
        ),
    ),
)
def case_data(request) -> CaseData:
    """
    A :class:`.CaseData` instance.

    """

    return request.param()
