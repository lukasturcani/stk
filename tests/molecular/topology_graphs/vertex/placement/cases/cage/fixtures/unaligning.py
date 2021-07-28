import pytest
import numpy as np
import stk

from ....case_data import CaseData


@pytest.fixture(
    params=(
        lambda: CaseData(
            vertex=stk.cage.UnaligningVertex(
                vertex=stk.cage.vertices._CageVertex(0, (1, 2, 3)),
            ),
            edges=(),
            building_block=stk.BuildingBlock(
                smiles='[Fe]',
                position_matrix=([0, 0, 0], ),
            ),
            position=np.array([1, 2, 3], dtype=np.float64),
            alignment_tests={},
            functional_group_edges={},
        ),
    ),
)
def unaligning(request) -> CaseData:
    return request.param()
