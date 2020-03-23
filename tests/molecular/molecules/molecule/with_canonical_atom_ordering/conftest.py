import pytest
import numpy as np
import stk

from .case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            molecule=stk.BuildingBlock(
                smiles='Br[C+2][N+]Cl',
                functional_groups=[stk.BromoFactory()],
            ),
            result=stk.BuildingBlock.init(
                atoms=(
                    stk.Cl(0),
                    stk.Br(1),
                    stk.C(2, 2),
                    stk.N(3, 1),
                ),
                bonds=(
                    stk.Bond(stk.Br(1), stk.C(2, 2), 1),
                    stk.Bond(stk.C(2, 2), stk.N(3, 1), 1),
                    stk.Bond(stk.N(3, 1), stk.Cl(0), 1),
                ),
                position_matrix=np.array([
                ]),
                functional_groups=(
                    stk.Bromo(
                        bromine=stk.Br(1),
                        atom=stk.C(2, 2),
                        bonders=(stk.C(2, 2), ),
                        deleters=(stk.Br(1), ),
                    ),
                ),
                placer_ids=(2, ),
            )
        ),
    ),
)
def case_data(request):
    return request.param
