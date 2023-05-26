import numpy as np
import pytest
import stk

from ....case_data import CaseData


@pytest.fixture(
    params=(
        lambda: CaseData(
            vertex=stk.cof.UnaligningVertex(
                id=0,
                position=(1, 2, 3),
            ),
            edges=(),
            building_block=stk.BuildingBlock(
                smiles="[Fe]",
                position_matrix=([0, 0, 0],),
            ),
            position=np.array([1, 2, 3], dtype=np.float64),
            alignment_tests={},
            functional_group_edges={},
        ),
    ),
)
def unaligning(request) -> CaseData:
    return request.param()
