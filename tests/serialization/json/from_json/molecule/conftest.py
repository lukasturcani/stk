import pytest
import numpy as np
import stk

from .case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            dejsonizer=stk.MoleculeDejsonizer(),
            json={
                'a': (
                    (35, 0),
                    (6, 2),
                    (6, 2),
                    (35, 0),
                ),
                'b': (
                    (0, 1, 1, (0, 0, 0)),
                    (1, 2, 1, (0, 0, 0)),
                    (2, 3, 1, (0, 0, 0)),
                ),
                'InChI': 'InChI=1S/C2Br2/c3-1-2-4/q+4',
                'InChIKey': 'UWAHASCVLDBPQQ-UHFFFAOYSA-N',
            },
            position_matrix=np.array([
                [0., 0., 0.],
                [1., 0., 0.],
                [2., 0., 0.],
                [3., 0., 0.],
            ]),
            molecule=stk.BuildingBlock(
                smiles='Br[C+2][C+2]Br',
            ).with_position_matrix(np.array([
                [0., 0., 0.],
                [1., 0., 0.],
                [2., 0., 0.],
                [3., 0., 0.],
            ])),
        ),
    ),
)
def case_data(request):
    """
    A :class:`.CaseData` instance.

    """

    return request.param
